#Packages needed
from scipy import io
import pandas as pd
import os
from tqdm import tqdm
import numpy as np
from pathlib import Path
from itertools import combinations
import random
import math
import matplotlib.pyplot as plt
from .general_utility_functions import check_input
random.seed(2022)

def count_doublets(df, output_prefix = None, dataset_name = None,
                   save_all_singlet_categories = False, save_plot_umi = False,
                   umi_cutoff_method = "ratio", umi_cutoff_ratio=3 / 4e5, umi_cutoff_percentile = None, min_umi_good_data_cutoff=2,
                   umi_diff_threshold=50, dominant_threshold=10):
    # TODO: write check_input function
    check_input(df, output_path = output_prefix, dataset_name = dataset_name,
                save_all_singlet_categories = save_all_singlet_categories,
                save_plot_umi = save_plot_umi,
                umi_cutoff_method = umi_cutoff_method)
    if(output_prefix == None):
        print("Note: None of the output files will be saved since no output path was provided.")
    else:
        output_prefix = os.path.join(output_prefix)
        if dataset_name:
            output_prefix = os.path.join(output_prefix, dataset_name + "_")

    df_sum = df['sample'].value_counts()
    all_samples = df['sample'].unique()
    print("INFO: Raw data counts: ")
    print(df_sum)


    # adjust UMI cutoff based on reads
    # drop cells that do not pass UMI cutoff

    #NEW Renaming good_data_ls to singlets_list
    singlets_list = []
    all_singlet_stats = []
    for cur_sample_num in df["sample"].unique():
        cur_singlet_stats = []

        # subset by sample
        cur_sample_df = df[df["sample"] == cur_sample_num]
        total_cells = len(cur_sample_df["cellID"].unique())
        print("Total cells for sample {}: {}".format(cur_sample_num, total_cells))

        #Counting nUMI for each barcode, cell ID combination
        cur_freq_df = cur_sample_df.groupby(['cellID', 'barcode', 'sample']).size().reset_index()
        cur_freq_df.columns = ['cellID', "barcode", 'sample', "nUMI"]
        
        if umi_cutoff_method == "percentile":
            print("INFO: Using percentile filtering by choosing barcodes whose UMI counts is in top {} percentile".format(umi_cutoff_percentile))
            calculated_umi_cutoff = math.ceil(np.percentile(cur_freq_df["nUMI"], umi_cutoff_percentile))
            cur_umi_adjusted_cutoff = max(calculated_umi_cutoff, min_umi_good_data_cutoff)
        
        elif umi_cutoff_method == "ratio":
            print("INFO: Using ratio based filtering.")
            calculated_umi_cutoff = round(umi_cutoff_ratio * cur_sample_df.shape[0])
            cur_umi_adjusted_cutoff = max(calculated_umi_cutoff, min_umi_good_data_cutoff)
        else:
            raise ValueError("The umi_cutoff_method chosen has to be percentile or ratio for RNAseq based barcodes")
        
        # Plotting the distribution of all data
        if(save_plot_umi):
            plt.hist(cur_freq_df["nUMI"], bins=30, edgecolor='black', alpha=0.7)

            # Add a red line for the 10th percentile
            plt.axvline(cur_umi_adjusted_cutoff, color='red', linestyle='dashed', linewidth=2)
            plt.title(f'UMI cutoff: {cur_umi_adjusted_cutoff} UMI')
            cur_umi_plot_file = output_prefix + "sample_" + str(cur_sample_num) + ".png"
            plt.savefig(cur_umi_plot_file)
        print("Current Sample Adjusted UMI cutoff: {}".format(cur_umi_adjusted_cutoff))

        #Filtering using UMI cutoff
        cur_good_data = cur_freq_df[cur_freq_df['nUMI'] >= cur_umi_adjusted_cutoff].reset_index(drop=True)
        singlets_list.append(cur_good_data)

        low_umi_cells_removed = total_cells - len(cur_good_data["cellID"].unique())
        
        # get fatemap barcode count for each cellID
        cellID2fatemap_dict, cellID2fatemap_count_dict = \
            generate_fatemap_barcode_counts_for_cellID(cur_good_data)

        # get barcodes present in multiple samples
        multilane_barcodes = get_multilane_barcodes(cur_good_data)

        # get first round of results
        singlets, multiplets = define_singlets_and_multiplets_based_on_barcode_counts(cellID2fatemap_count_dict,
                                                                                        cur_good_data,
                                                                                        multilane_barcodes)
        cur_singlet_stats.append(len(singlets))

        # salvage false multiplets
        barcode_dict, barcode_count_dict = generate_barcode_combo(multiplets, cur_good_data)
        two_barcode_singlets = extract_two_barcode_singlets(barcode_dict,
                                                              barcode_count_dict)
        cur_singlet_stats.append(len(two_barcode_singlets))

        # update singlets and multiplets
        singlets_step2 = list(set(singlets).union(set(two_barcode_singlets)))
        multiplets_step2 = list(set(multiplets).difference(set(two_barcode_singlets)))

        # recover singlets with dominant UMI
        UMI_thres_singlets = identify_singlets_with_dominant_UMI(multiplets_step2,
                                                                 cur_good_data,
                                                                 umi_diff_threshold,
                                                                 dominant_threshold)
        cur_singlet_stats.append(len(UMI_thres_singlets))
        singlets_step3 = list(set(singlets_step2).union(set(UMI_thres_singlets)))
        multiplets_step3 = list(set(multiplets_step2).difference(set(UMI_thres_singlets)))

        cur_singlet_stats.append(len(singlets_step3))

        cur_singlet_stats.append(len(multiplets_step3))
        print("Total Singlets: {}\nTotal Multiplets: {}".format(len(singlets_step3), len(multiplets_step3)))

        singlets_stat_df = pd.DataFrame(cur_singlet_stats).T
        singlets_stat_df.columns = ["single_sample_barcode_singlets", "multi_barcode_singlets",
                                    "dominant_umi_barcode_singlets", "total_singlets", "total_undetermined"]
        singlets_stat_df["low_umi_cells_removed"] = low_umi_cells_removed
        singlets_stat_df["total_cells"] = total_cells
        if(dataset_name):
            singlets_stat_df["dataset"] = dataset_name
        singlets_stat_df["sample"] = cur_sample_num

        #Keeping track of stats for all samples
        all_singlet_stats.append(singlets_stat_df)
        cur_good_data["label"] = cur_good_data["cellID"].apply(
            lambda x: "Singlet" if x in singlets_step3 else "Multiplet")


        #This to create pairs of singlets for simulating doublets - so may remove

        # generate singlet pair for each sample for doublet simulation
        # cur_singlets = cur_good_data[cur_good_data["label"] == "Singlet"]["cellID"].values
        # singlets_combo_for_doublets_simulation = [", ".join(map(str, comb)) for comb in combinations(cur_singlets, 2)]
        # random_index = random.sample(range(0, len(singlets_combo_for_doublets_simulation)), len(cur_singlets))
        # cur_out_file = output_prefix + "_{}_singlet_pairs.csv".format(cur_sample_num)
        # with open(cur_out_file, "w+") as fp:
        #     fp.writelines("%s\n" % singlets_combo_for_doublets_simulation[idx] for idx in random_index)

        # with open(simulated_doublets_file, "w+") as fp:
        #     singlets_combo_for_doublets_simulation = [", ".join(map(str, comb)) for comb in combinations(singlets_step3, 2)]
        #     random_index = random.sample(range(0, len(singlets_combo_for_doublets_simulation)), len(singlets_step3))
        #     fp.writelines("%s\n" % singlets_combo_for_doublets_simulation[idx] for idx in random_index)
        # initialize output files
        if(save_all_singlet_categories):
            singlets_file = output_prefix + "_singlets_all.txt"
            singlets_single_barcode_file = output_prefix + f"_{cur_sample_num}_" + "_single_barcode_singlets.txt"
            singlets_multi_barcode_file = output_prefix + f"_{cur_sample_num}_" + "_multi_barcode_singlets.txt"
            singlets_dominant_umi_file = output_prefix + f"_{cur_sample_num}_" + "_dominant_umi_singlets.txt"
            multiplets_file = output_prefix + f"_{cur_sample_num}_" + "_multiplets.txt"
            with open(singlets_single_barcode_file, "w+") as fp:
                fp.writelines("%s\n" % l for l in singlets)
            with open(singlets_multi_barcode_file, "w+") as fp:
                fp.writelines("%s\n" % l for l in two_barcode_singlets)
            with open(singlets_dominant_umi_file, "w+") as fp:
                fp.writelines("%s\n" % l for l in UMI_thres_singlets)
            with open(singlets_file, "w+") as fp:
                fp.writelines("%s\n" % l for l in singlets_step3)
            with open(multiplets_file, "w+") as fp:
                fp.writelines("%s\n" % l for l in multiplets_step3)
            singlet_stats_file = output_prefix + f"_{cur_sample_num}_" + "singlets_stats.csv"
            singlets_stat_df.to_csv(singlet_stats_file, index=False)
        
    good_data = pd.concat(singlets_list, ignore_index=True)
    all_stats = pd.concat(all_singlet_stats, ignore_index=True)


    #TODO Integrate multi-sample singlets into above - singlets_step4
    # concat all good data to identify multi sample multi fatemap barcode
    # This is a very rare case for singlets to be specific to be cross sample multibarcode
    # They are calculated but not assigned back to each sample, but is potentially an improvement in the future
    
    good_data = good_data.sort_values(by=["sample", 'cellID', 'barcode', 'nUMI']).reset_index(drop=True)
    multilane_barcodes = get_multilane_barcodes(good_data)
    if(output_prefix):
        all_cross_sample_cell_id_file = output_prefix + "_all_sample_two_barcode_cell_id.txt"
        with open(all_cross_sample_cell_id_file, "w+") as fp:
            fp.writelines("%s\n" % l for l in multilane_barcodes)
    # barcode_dict, barcode_count_dict = generate_barcode_combo(multilane_barcodes,
    #                                                                 good_data)
    # all_sample_two_barcode_singlets = extract_two_barcode_singlets(barcode_dict,
    #                                                                  barcode_count_dict)
    # all_sample_two_barcode_singlets_barcode_file = output_prefix + "_all_sample_two_barcode_singlets.txt"
    # with open(all_sample_two_barcode_singlets_barcode_file, "w+") as fp:
    #     fp.writelines("%s\n" % l for l in all_sample_two_barcode_singlets)

    return good_data,  all_stats

#Function to identify singlets with more than one barcode but one of them is dominant
def identify_singlets_with_dominant_UMI(multiplets_step2,
                                        good_data,
                                        umi_diff_threshold,
                                        dominant_threshold):
    UMI_thres_singlets = []
    for cur_multiplet in multiplets_step2:
        cur_df = good_data[good_data["cellID"] == cur_multiplet]
        cur_UMI_counts = cur_df['nUMI'].values
        cur_median_UMI = np.median(cur_UMI_counts)
        if (np.max(cur_UMI_counts) - cur_median_UMI) > umi_diff_threshold:
            # check if there is more than 1 dominant
            dominant_count = 0
            for cur_UMI in cur_UMI_counts:
                cur_diff = cur_UMI - cur_median_UMI
                if cur_diff > umi_diff_threshold or cur_UMI > dominant_threshold:
                    dominant_count += 1
            if dominant_count == 1:
                UMI_thres_singlets.append(cur_multiplet)

    return UMI_thres_singlets


def extract_two_barcode_singlets(barcode_dict, barcode_count_dict):
    two_barcode_singlets_count = 0
    two_barcode_singlets = []
    for key, value in barcode_count_dict.items():
        if value >= 2:
            two_barcode_singlets_count += value
            two_barcode_singlets.extend(barcode_dict[key])
            if barcode_dict[key] in two_barcode_singlets:
                print(barcode_dict[key])
    print("All singlets identified with multiple barcodes are unique? {}".format(check_all_unique(two_barcode_singlets)))
    # print("Count of salvaged cells: {}\nProportion of salvaged cells: {}".format(two_barcode_singlets_count,
    # two_barcode_singlets_count/len(good_data["cellID"].unique())))
    return two_barcode_singlets

def define_singlets_and_multiplets_based_on_barcode_counts(cellID2fatemap_count_dict,
                                                             good_data,
                                                             blacklist_barcodes):
    singlets = []
    multiplets = []
    for sample_id, cur_dict in cellID2fatemap_count_dict.items():
        cur_sample_count = 0
        for barcode, count in cur_dict.items():
            # if the barcode is present in more than one lane, remove it
            if barcode in blacklist_barcodes:
                continue
            if count == 1:
                cur_sample_count += 1
                singlets.append(barcode)
            else:
                multiplets.append(barcode)
    return singlets, multiplets

#Function to get singlets in multiple samples
def get_multilane_barcodes(good_data):
    multilane_barcodes = []
    for cur_barcode in good_data["cellID"].unique():
        cur_df = good_data[good_data["cellID"] == cur_barcode]
        if (len(cur_df["sample"].unique()) > 1):
            multilane_barcodes.append(cur_barcode)
    return multilane_barcodes

# Check how may of these are from the same sample and how many are from different samples
def generate_barcode_combo(multiplets, good_data):
    barcode_dict = {}
    barcode_count_dict = {}
    # upon close examination, this should account for both multibarcode singlet within sample and cross samples
    # this is because we are counting the unique combo of barcode in good_data regardless of samples
    for cur_multiplet in multiplets:
        cur_df = good_data[good_data["cellID"] == cur_multiplet]
        # only consider when more than 1 fate barcode appears
        n_unique_fatemap_barcodes = len(cur_df["barcode"].unique())

        # a bit redundant since multiplets by definition should have more than one unique fatemap barcode associated
        if n_unique_fatemap_barcodes < 2:
            continue
        # create unique key for each fatemap barcode ID
        barcode_combo = "_".join(sorted(cur_df["barcode"].unique()))
        if barcode_combo not in barcode_dict:
            barcode_dict[barcode_combo] = [cur_multiplet]
            barcode_count_dict[barcode_combo] = 1
        else:
            barcode_dict[barcode_combo].append(cur_multiplet)
            barcode_count_dict[barcode_combo] += 1
    return barcode_dict, barcode_count_dict

#Won't be needed for singletCode - make separate function
def generate_fatemap_barcode_counts_for_cellID(good_data):
    all_samples = good_data['sample'].unique()
    cellID2fatemap_dict = {}
    cellID2fatemap_count_dict = {}

    # initialize the count dictionaries
    for i in all_samples:
        cellID2fatemap_dict[i] = {}
        cellID2fatemap_count_dict[i] = {}

    # loop through each row
    for index, row in tqdm(good_data.iterrows(), total=good_data.shape[0]):
        cur_cellID = row['cellID']
        cur_fateID = row['barcode']
        cur_sample_num = row['sample']

        # only look at samples 1 and 2 for now
        # if(cur_sample_num not in all_samples):
        #     continue
        if cur_cellID not in cellID2fatemap_dict[cur_sample_num].keys():
            cellID2fatemap_dict[cur_sample_num][cur_cellID] = [cur_fateID]
            cellID2fatemap_count_dict[cur_sample_num][cur_cellID] = 1
        else:
            # only add if the fateID is not present in the dict
            if cur_fateID not in cellID2fatemap_dict[cur_sample_num][cur_cellID]:
                cellID2fatemap_dict[cur_sample_num][cur_cellID].append(cur_fateID)
                cellID2fatemap_count_dict[cur_sample_num][cur_cellID] += 1
    return cellID2fatemap_dict, cellID2fatemap_count_dict


def check_all_unique(target):
    return len(list(set(target))) == len(target)
