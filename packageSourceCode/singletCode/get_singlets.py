#Module which has the the main functions
from .count_doublets_utils_copy import count_doublets

def get_singlets(sample_sheet, output_path = None, dataset_name = None, sample_type = "RNA",
                save_all_singlet_categories = False, save_plot_umi = False, 
                umi_cutoff_method = "ratio", umi_cutoff_ratio = 3/4e5, 
                umi_cutoff_percentile= None,  min_umi_cutoff =2, umi_diff_threshold = 50,
                umi_dominant_threshold = 10):
    """
    Function that inputs the sample sheet and other parameters and runs it through count_doublets to get a list of singlets in the sample.
    :param: sample_sheet: a dataframe that contains 3 columns: cellID, barcode, sample and if a row is repeated, it is assumed to reflect the UMI associated with a barcode in this cell (identified by the cellID). 
    :dataset_name: This is the name of the dataset being analysed. It will be in the name of all saved files and be a column in the singlet_stats sheet returned
    :param: output_path: the path to store any output files, including plots to show UMI distribution and what the umi_cutoff used is, csv files containing singlets of different categories. If None, then the list of singlets will be returned but it won't contain information about what category of singlet each cell is.
    :param: sample_type: specify if the barcodes are from RNAseq or ATACseq data. If it is from ATACseq data, any barcode with more than 2 UMI will be filtered out since there can be atmost 2 copies of the barcode in the genome and hence, atmost 2 UMI associated with it. Also, please ensure the input_sheet contains UMI information, and not reads.  
    :param: save_all_singlet_categories: If true, then singlets of each category are saved separately in csv files along with all singlets and all non-singlets - a total of 5 files - at the path inputted in output_prefix
    :param: save_plot_umi: If true, then plots showing UMI distribution indicating the UMI cutoff used will be saved for each sample in the path inputted in output_prefix.
    :param: umi_cutoff_method: Specify if quality control for barcodes using UMI counts should be based on "ratio" or "percentile". If it is ratio, then the umi_cutoff will be calculated by multiplying umi_cutoff_ratio (default: 3/4e5) with the total UMI count for sample. If it is percentile, the umi_cutoff_percentile will be used to determine the umi_cutoff.
    :param: umi_cutoff_ratio: this is the ratio used to determine the umi_cutoff if umi_cutoff_method is "ratio"
    :param: umi_cutoff_percentile: If umi_cutoff_method is "percentile", then the umi_cutoff will be the minimum UMI count required to be in top umi_cutoff_percentile'th percentile. 
    :param: min_umi_cutoff: This is the absolute minimum number of UMIs that need to be associated with a barcode for it to be considered a barcode. However, the actual umi_cutoff used will be the greater of min_umi_cutoff and the cutoff calculated using umi_cutoff_method.    
    :param: umi_dominant_threshold: The minimum UMI count to be associated with a barcode for it to be considered to be a potential dominant barcode in a cell.
    :param: umi_diff_threshold: This is the minimum difference between UMI counts associated with a potential dominant barcode present within a cell and the median UMI count of all barcodes associated with the cell. If a cell has only one dominant barcode, it will be counted as a singlet.

    """
    singlet_list, singlet_stats = count_doublets(
            sample_sheet,
            output_prefix=output_path,
            save_all_singlet_categories = save_all_singlet_categories,
            save_plot_umi = save_plot_umi,
            umi_cutoff_ratio=umi_cutoff_ratio,
            umi_cutoff_percentile = umi_cutoff_percentile,
            umi_diff_threshold=umi_diff_threshold,
            dominant_threshold=umi_dominant_threshold,
            min_umi_good_data_cutoff=min_umi_cutoff,
            sample_type = sample_type, 
            umi_cutoff_method = umi_cutoff_method,
            dataset_name = dataset_name
    )
    
    return singlet_list, singlet_stats

