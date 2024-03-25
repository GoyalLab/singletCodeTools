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

random.seed(2022)


def count_doublets(input_file, output_prefix, umi_cutoff_ratio=3 / 4e5, umi_diff_threshold=50, dominant_threshold=10,
                   overwrite=True, min_umi_good_data_cutoff=2):
    ext = os.path.splitext(input_file)[1].lower()
    if ext == ".csv":
        df = pd.read_csv(input_file)
    elif ext == ".tsv" or ext == ".txt":
        df = pd.read_csv(input_file, sep='\t')
    else:
        print("Does not recognize extension. Quiting...")
        return None

    # TODO: write a merge function at the end
    singlets_file = output_prefix + "_singlets_all.txt"

    if (Path(singlets_file).is_file() or Path(singlets_file).is_file()) and overwrite == False:
        print("Output files exist for input: {}\nUse `overwrite` parameter to overwrite".format(input_file))
        return None

    # some sampleNum is written as SampleNum
    # force conversion
    if df.shape[1] == 3:
        df.columns = ['cellID', 'BC50StarcodeD8', 'sampleNum']
    # some data has four columns, one of which is UMI
    elif df.shape[1] == 4:
        df = df.drop(columns=['UMI'])
        df.columns = ['cellID', 'BC50StarcodeD8', 'sampleNum']

    # check data
    else:
        print("More than 3 columns. Quiting...")
        return None
    # subset for cellTag
    if "cellTag" in input_file:
        # strings_to_keep = ['d2-RNA-5', 'd2-ATAC-5', 'B4D21-RNA-r2-1', 'B4D21-ATAC-r2-1']
        strings_to_keep = ['d2-RNA-5', 'B4D21-RNA-r2-1']
        df = df[df['sampleNum'].isin(strings_to_keep)]

    df_sum = df['sampleNum'].value_counts()
    all_samples = df['sampleNum'].unique()
    print("INFO: Raw data counts")
    print(df_sum)


    # adjust UMI cutoff based on reads
    # drop cells that do not pass UMI cutoff
    good_data_ls = []
    for cur_sample_num in df["sampleNum"].unique():
        cur_singlet_stats = []

        # subset by sampleNum
        cur_sample_df = df[df["sampleNum"] == cur_sample_num]
        total_cells = len(cur_sample_df["cellID"].unique())

        cur_freq_df = cur_sample_df.groupby(['cellID', 'BC50StarcodeD8', 'sampleNum']).size().reset_index()
        cur_freq_df.columns = ['cellID', "fatemapID", 'sampleNum', "nUMI"]

        # if it is atac data, then the sample id may have ATAC in it, for instance for the cellTag dataset
        # or it is the watermelon dataset
        # in this case we remove all cells with greater than 2 UMI
        # if ("watermelon" in input_file) or ("ATAC" in cur_sample_num):
        if "ATAC" in str(cur_sample_num):
            print(cur_sample_num)
            cur_freq_df = cur_freq_df[cur_freq_df['nUMI'] <= 2].reset_index(drop=True)
            calculated_umi_cutoff = 1
            cur_umi_adjusted_cutoff = 1
        elif umi_cutoff_ratio >= 1:
            cur_umi_plot_file = output_prefix + "_sample_" + str(cur_sample_num) + ".png"
            print("INFO: Using percentile filtering. e.g. Pass in 1 for filtering out the cells with lowest 1 percent "
                  "UMI count.")
            calculated_umi_cutoff = math.ceil(np.percentile(cur_freq_df["nUMI"], umi_cutoff_ratio))
            # Plotting the distribution of all data
            plt.hist(cur_freq_df["nUMI"], bins=30, edgecolor='black', alpha=0.7)

            # Add a red line for the 10th percentile
            plt.axvline(calculated_umi_cutoff, color='red', linestyle='dashed', linewidth=2)
            plt.title(f'{umi_cutoff_ratio}th percentile: {calculated_umi_cutoff} UMI')
            plt.savefig(cur_umi_plot_file)
            cur_umi_adjusted_cutoff = max(calculated_umi_cutoff, min_umi_good_data_cutoff)
        else:
            print("INFO: Using raio based filtering.")
            calculated_umi_cutoff = round(umi_cutoff_ratio * cur_sample_df.shape[0])
            cur_umi_adjusted_cutoff = max(calculated_umi_cutoff, min_umi_good_data_cutoff)
        print("Current Sample Adjusted UMI cutoff: {}".format(cur_umi_adjusted_cutoff))

        cur_good_data = cur_freq_df[cur_freq_df['nUMI'] >= cur_umi_adjusted_cutoff].reset_index(drop=True)
        good_data_ls.append(cur_good_data)
        low_umi_cells_removed = total_cells - len(cur_good_data["cellID"].unique())
        # get fatemap barcode count for each cellID
        cellID2fatemap_dict, cellID2fatemap_count_dict = \
            generate_fatemap_barcode_counts_for_cellID(cur_good_data)

        # get multilane barcodes
        multilane_barcodes = get_multilane_barcodes(cur_good_data)

        # get first round of results
        singlets, multiplets = define_singlets_and_multiplets_based_on_fatemapID_counts(cellID2fatemap_count_dict,
                                                                                        cur_good_data,
                                                                                        multilane_barcodes)
        cur_singlet_stats.append(len(singlets))

        # salvage false multiplets
        fatemapID_dict, fatemapID_count_dict = generate_fatemapID_combo(multiplets, cur_good_data)

        two_barcode_singlets = extract_two_fatemapID_singlets(fatemapID_dict,
                                                              fatemapID_count_dict)
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
        singlets_stat_df["dataset"] = input_file.split("/")[-3]
        singlet_stats_file = output_prefix + f"_{cur_sample_num}_" + "singlets_stats.csv"
        singlets_stat_df.to_csv(singlet_stats_file, index=False)
        cur_good_data["label"] = cur_good_data["cellID"].apply(
            lambda x: "Singlet" if x in singlets_step3 else "Multiplet")

        # generate singlet pair for each sample for doublet simulation
        cur_singlets = cur_good_data[cur_good_data["label"] == "Singlet"]["cellID"].values
        singlets_combo_for_doublets_simulation = [", ".join(map(str, comb)) for comb in combinations(cur_singlets, 2)]
        random_index = random.sample(range(0, len(singlets_combo_for_doublets_simulation)), len(cur_singlets))
        cur_out_file = output_prefix + "_{}_singlet_pairs.csv".format(cur_sample_num)
        with open(cur_out_file, "w+") as fp:
            fp.writelines("%s\n" % singlets_combo_for_doublets_simulation[idx] for idx in random_index)

        # with open(simulated_doublets_file, "w+") as fp:
        #     singlets_combo_for_doublets_simulation = [", ".join(map(str, comb)) for comb in combinations(singlets_step3, 2)]
        #     random_index = random.sample(range(0, len(singlets_combo_for_doublets_simulation)), len(singlets_step3))
        #     fp.writelines("%s\n" % singlets_combo_for_doublets_simulation[idx] for idx in random_index)
        # initialize output files
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

    # concat all good data to identify multi sample multi fatemap barcode
    # This is a very rare case for singlets to be specific to be cross sample multibarcode
    # They are calculated but not assigned back to each sample, but is potentially an improvement in the future
    good_data = pd.concat(good_data_ls)
    good_data = good_data.sort_values(by=["sampleNum", 'cellID', 'fatemapID', 'nUMI']).reset_index(drop=True)
    multilane_barcodes = get_multilane_barcodes(good_data)
    all_cross_sample_cell_id_file = output_prefix + "_all_sample_two_barcode_cell_id.txt"
    with open(all_cross_sample_cell_id_file, "w+") as fp:
        fp.writelines("%s\n" % l for l in multilane_barcodes)
    fatemapID_dict, fatemapID_count_dict = generate_fatemapID_combo(multilane_barcodes,
                                                                    good_data)
    all_sample_two_barcode_singlets = extract_two_fatemapID_singlets(fatemapID_dict,
                                                                     fatemapID_count_dict)
    all_sample_two_barcode_singlets_barcode_file = output_prefix + "_all_sample_two_barcode_singlets.txt"
    with open(all_sample_two_barcode_singlets_barcode_file, "w+") as fp:
        fp.writelines("%s\n" % l for l in all_sample_two_barcode_singlets)
# [1,2,2,3,60,190]
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


def extract_two_fatemapID_singlets(fatemapID_dict, fatemapID_count_dict):
    two_barcode_singlets_count = 0
    two_barcode_singlets = []
    for key, value in fatemapID_count_dict.items():
        if value >= 2:
            two_barcode_singlets_count += value
            two_barcode_singlets.extend(fatemapID_dict[key])
            if fatemapID_dict[key] in two_barcode_singlets:
                print(fatemapID_dict[key])
    print("All singlets identified are unique? {}".format(check_all_unique(two_barcode_singlets)))
    # print("Count of salvaged cells: {}\nProportion of salvaged cells: {}".format(two_barcode_singlets_count,
    # two_barcode_singlets_count/len(good_data["cellID"].unique())))
    return two_barcode_singlets


# TODO
# Check how may of these are from the same sample and how many are from different samples
def generate_fatemapID_combo(multiplets, good_data):
    fatemapID_dict = {}
    fatemapID_count_dict = {}
    # upon close examination, this should account for both multibarcode singlet within sample and cross samples
    # this is because we are counting the unique combo of fatemapID in good_data regardless of samples
    for cur_multiplet in multiplets:
        cur_df = good_data[good_data["cellID"] == cur_multiplet]
        # only consider when more than 1 fate barcode appears
        n_unique_fatemap_barcodes = len(cur_df["fatemapID"].unique())

        # a bit redundant since multiplets by definition should have more than one unique fatemap barcode associated
        if n_unique_fatemap_barcodes < 2:
            continue
        # create unique key for each fatemap barcode ID
        fatemapID_combo = "_".join(sorted(cur_df["fatemapID"].unique()))
        if fatemapID_combo not in fatemapID_dict:
            fatemapID_dict[fatemapID_combo] = [cur_multiplet]
            fatemapID_count_dict[fatemapID_combo] = 1
        else:
            fatemapID_dict[fatemapID_combo].append(cur_multiplet)
            fatemapID_count_dict[fatemapID_combo] += 1
    return fatemapID_dict, fatemapID_count_dict


def define_singlets_and_multiplets_based_on_fatemapID_counts(cellID2fatemap_count_dict,
                                                             good_data,
                                                             blacklist_barcodes):
    singlets = []
    multiplets = []
    print("Total cells: {}".format(len(good_data["cellID"])))
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
        print("Sample {} singlet: {}".format(sample_id, cur_sample_count))
    print("Total Singlets: {}\nTotal Multiplets: {}".format(len(singlets), len(multiplets)))
    return singlets, multiplets


def get_multilane_barcodes(good_data):
    multilane_barcodes = []
    for cur_barcode in good_data["cellID"].unique():
        cur_df = good_data[good_data["cellID"] == cur_barcode]
        if (len(cur_df["sampleNum"].unique()) > 1):
            multilane_barcodes.append(cur_barcode)
    return multilane_barcodes


def generate_fatemap_barcode_counts_for_cellID(good_data):
    all_samples = good_data['sampleNum'].unique()
    cellID2fatemap_dict = {}
    cellID2fatemap_count_dict = {}

    # initialize the count dictionaries
    for i in all_samples:
        cellID2fatemap_dict[i] = {}
        cellID2fatemap_count_dict[i] = {}

    # loop through each row
    for index, row in tqdm(good_data.iterrows(), total=good_data.shape[0]):
        cur_cellID = row['cellID']
        cur_fateID = row['fatemapID']
        cur_sample_num = row['sampleNum']

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
