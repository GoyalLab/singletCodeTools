#Module to define the general utility functions.
import pandas as pd
import os
def check_sample_sheet(data_frame):
    '''
    Function to check if the dataframe can be used as input to get_singlets function. It checks if the three columns - cellID, barcode, sample needed are present and in same order. If this dataframe can be used for get_singlets function, then a statement will be printed to confirm that
    
    Args: 
        sample_sheet: A dataframe that contains your sample sheet.
    
    Returns:
        None
    '''
    column_names = data_frame.columns
    column_names = column_names.tolist()
    required_columns = ["cellID","barcode", "sample"]
    if column_names != required_columns:
        #Test
        raise ValueError("The column names of the sample sheet must be cellID, barcode and sample, in that order.")
    print("The sample sheet provided can be used as input to get_singlets to get a list of singlets identified.")
    return

def check_input(sample_sheet, output_path = None, dataset_name = None, 
                save_all_singlet_categories = None, save_plot_umi = None,
                umi_cutoff_method = None):
                

    if(output_path):
        flag = os.path.exists(output_path)
        if(not flag):
            raise ValueError("Verify that the provided output_path exists.")
    #checking sample sheet
    column_names = sample_sheet.columns
    column_names= column_names.tolist()
    required_columns = ["cellID", "barcode", "sample"]
    if column_names != required_columns:
        raise ValueError("The column names of the sample sheet must be cellID, barcode and sample, in that order.")
    
    #checking if output_prefix is provided if save is true
    if (save_plot_umi or save_all_singlet_categories):
        if(not output_path):
            raise ValueError("The option to save files has been set to true, but no output_path has been provided.")
        
    if (umi_cutoff_method != "ratio" and umi_cutoff_method != "percentile"):
        raise ValueError("Please choose either ratio or percentile method to determine UMI cutoff for you RNAseq dataset.")
    return