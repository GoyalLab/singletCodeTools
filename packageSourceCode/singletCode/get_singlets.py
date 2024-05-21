#Module which has the the main functions
from .count_doublets_utils_copy import count_doublets

def get_singlets(sample_sheet, output_path = None, dataset_name = None,
                save_all_singlet_categories = False, save_plot_umi = False, 
                umi_cutoff_method = "ratio", umi_cutoff_ratio = 3/4e5, 
                umi_cutoff_percentile= None,  min_umi_cutoff =2, umi_diff_threshold = 50,
                umi_dominant_threshold = 10):
    """
    Function that inputs the sample sheet and other parameters and runs it through count_doublets to get a list of 
    singlets in the sample. If a row is repeated, it is assumed to reflect the UMI associated with a barcode in 
    this cell (identified by the cellID).

    Args:
        sample_sheet (DataFrame): A dataframe that contains 3 columns: cellID, barcode, sample.
        dataset_name (str): The name of the dataset being analysed. It will be in the name of all saved files 
            and be a column in the singlet_stats sheet returned.
        output_path (str, optional): The path to store any output files, including plots to show UMI distribution 
            and what the umi_cutoff used is, csv files containing singlets of different categories. If None, 
            then the list of singlets will be returned but it won’t contain information about what category of 
            singlet each cell is. Defaults to None.
        save_all_singlet_categories (bool, optional): If true, then singlets of each category are saved separately 
            in csv files along with all singlets and all non-singlets. Defaults to False.
        save_plot_umi (bool, optional): If true, then plots showing UMI distribution indicating the UMI cutoff used 
            will be saved for each sample. Defaults to False.
        umi_cutoff_method (str, optional): Specify if quality control for barcodes using UMI counts should be based 
            on "ratio" or "percentile". Defaults to 'ratio'.
        umi_cutoff_ratio (float, optional): The ratio used to determine the umi_cutoff if umi_cutoff_method is 
            "ratio". Defaults to 3/4e5.
        umi_cutoff_percentile (float, optional): If umi_cutoff_method is "percentile", then the umi_cutoff will be 
            the minimum UMI count required to be in the top umi_cutoff_percentile'th percentile. There is no default and if umi_cutoff_method is set to "percentile", then manually set this parameter.
        min_umi_cutoff (int, optional): This is the absolute minimum number of UMIs that need to be associated with 
            a barcode for it to be considered a barcode. However, the actual umi_cutoff used will be the greater 
            of min_umi_cutoff and the cutoff calculated using umi_cutoff_method. Defaults to 10.
        umi_dominant_threshold (int, optional): The minimum UMI count to be associated with a barcode for it to be 
            considered to be a potential dominant barcode in a cell. Defaults to 5.
        umi_diff_threshold (int, optional): This is the minimum difference between UMI counts associated with a 
            potential dominant barcode present within a cell and the median UMI count of all barcodes associated 
            with the cell. If a cell has only one dominant barcode, it will be counted. Defaults to 2.

    Returns:
        tuple: A 2-tuple containing:
            - pandas.DataFrame: A dataframe which contains all unique cell ID and barcode combinations in the data along with singlet assignment to each cell ID.
            - pandas.DataFrame: A dataframe which contains the statistics for total singlets, different categories of singlets, and cells removed due to low UMI counts.
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
            umi_cutoff_method = umi_cutoff_method,
            dataset_name = dataset_name
    )
    
    return singlet_list, singlet_stats

