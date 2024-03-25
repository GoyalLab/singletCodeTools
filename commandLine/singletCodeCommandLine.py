"""
SingletCode helpful scripts

Description:
    This script contains two modules.
    The count module generates the singlet files for the input data sheet.
    The watermelon module uses the miseq dial out files to create the cell ID, 
    barcode and sample file and can be used as input to singlet code module

Usage:
    python singletCode.py count -i /path/to/input.txt -o /path/to/output
    OR
    python file_processor.py --input_file /path/to/input.txt --output_file /path/to/output.txt
    OR python3 singletCode.py watermelon -i /path/to/fastq/files -o path/to/save/csv/file 
    -s path/sample/sheet -use10X False -input10X path/to/barcodes/tsv 

Options for count module:
    -i --input_file    Specify the path to the input barcode file.
    -o --output_prefix   Specify the path to the output prefix.
    -f --force  Force overwrite if output file already exists
    -u --cutoff UMI cutoff ratio

Options for watermelon module:
    -i --inputFolder Specify the path to the folder containing the fastq folders output from miseq
    -s --sampleSheet Specify the path to the sample sheet in .csv format that contains sample name and sample number such that it matches the names of the fastq files (expected naming format which is typical for miseq output files: sampleName_sampleNumber_L001_ReadNumber_001.fastq.gz; an example: Sample-1_S1_L001_R1_001.fastq.gz)
    -o --outputFolder Specify the path to save the output csv file containing the barcode, cell ID information
    -outputName Specify the name of output csv file
    --use10X Specify if a 10X object is provided which has the same cells as the ones in fastq; if provided, the cells in the fastq file will be filtered out if not prsent in the 10X object
    --input10X Path to the barcodes.tsv.gz or barcodes.tsv which contains cell IDs

Author:
    Ziyang Zhang (Charles)
    Keerthana M Arun
"""

from scipy import io
import pandas as pd
import os
import numpy as np
from pathlib import Path
import pathlib
import argparse
import shutil
import os
import numpy as np
import argparse
import shutil
from count_doublets_utils import *
from watermelonUtilityFunctions import checkInputs, processFastqFiles

################################
#                              #
#       Input Parameters       #
#                              #
################################

# data_root = "/projects/p31666/zzhang/doublet-bchmk/data/fatemap_data"
# keyword = "stepFourStarcodeShavedReads50"
# pattern = '**/{}*'.format(keyword)
# overwrite = True
# # output_prefix = None
# output_prefix = "/projects/p31666/zzhang/doublet-bchmk/data/couting_threshold/scale_factor_0.00003"
# # umi_cutoff_ratio = 20
# umi_cutoff_ratio = 3 / 4e5


if __name__ == "__main__":
    # Create the top-level parser
    parser = argparse.ArgumentParser(description='SingletCode helper scripts')
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # Create the parser for the "copy" command
    parser_count = subparsers.add_parser('count', help='Counting module')
    parser_count.add_argument('-i', '--input_file', type=str, help='Input barcode file', required=True)
    parser_count.add_argument('-o', '--out_prefix', type=str, help='Output prefix', required=True)
    parser_count.add_argument('-u', '--umi_cutoff_ratio', type=float, help='Cutoff ratio for UMI filter',
                              required=False, default=3 / 4e5)
    parser_count.add_argument('-d', '--umi_diff_threshold', type=int, help='Minimum difference in umi', required=False,
                              default=50)
    parser_count.add_argument('-m', '--dominant_threshold', type=int, help='Minimum dominant umi', required=False,
                              default=10)
    parser_count.add_argument('-g', '--min_umi_good_data_cutoff', type=int, help='Minimum umi filter threshold',
                              required=False, default=2)

    # Create the parser for the "rename" command
    # parser_rename = subparsers.add_parser('rename', help='Rename a file')
    # parser_rename.add_argument('-i', '--input_file', type=str, help='Input barcode file', required=True)
    # parser_rename.add_argument('-n', '--new_name', type=str, help='New file name', required=True)
    
    #Create parser for watermelon module
    #Creating the parser for processing watermelon fastq data

    parser_watermelon = subparsers.add_parser('watermelon', help='watermelon module')
    parser_watermelon.add_argument('-i', '--inputFolder', type=str,
                                    help="Path to folder containing fastq files", required=True)
    parser_watermelon.add_argument('-o', '--outputFolder', type=str,
                                   help='Output Path', required=True)
    parser_watermelon.add_argument('-s', '--sampleSheet', type=str,
                                   help='Path to csv file containing sample details', 
                                   required=True)
    parser_watermelon.add_argument('--outputName', default="watermelonBarcodeUMIcount.csv",
                                   type=str,
                                   help='Name of output csv file', required=False)
    parser_watermelon.add_argument('--use10X', type=bool,
                                   help='if 10X cell IDs to be used to match cells in fastq files', required=False, default=False)
    parser_watermelon.add_argument('--input10X', type=str,
                                   help='Path to cell IDs of 10X object', required=False)
                                   
    # Parse the arguments
    args = parser.parse_args()
    print("Arguments received:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")
    # Call the appropriate function based on the sub-command
    if args.command == 'count':
        input_file = args.input_file
        output_prefix = args.out_prefix
        umi_cutoff_ratio = args.umi_cutoff_ratio
        umi_diff_threshold = args.umi_diff_threshold
        dominant_threshold = args.dominant_threshold
        min_umi_good_data_cutoff = args.min_umi_good_data_cutoff
        count_doublets(
            input_file=input_file,
            output_prefix=output_prefix,
            umi_cutoff_ratio=umi_cutoff_ratio,
            umi_diff_threshold=umi_diff_threshold,
            dominant_threshold=dominant_threshold,
            min_umi_good_data_cutoff=min_umi_good_data_cutoff
        )
    elif args.command == 'trim':
        pass
    elif args.command == "watermelon":
        inputFolder = args.inputFolder
        outputFolder = args.outputFolder
        sampleSheet = args.sampleSheet
        outputFileName = args.outputName
        use10Xcells = args.use10X
        if args.use10X and not args.input10X:
            parser.error("use10X to is set to True, but path to barcodes.tsv (or barcodes.tsv.gz) was not provided.")
        elif not args.use10X and not args.input10X:
            checkInputs(inputFolder, outputFolder, sampleSheet)
            processFastqFiles(inputFolder, outputFolder, sampleSheet, outputFileName)
        else:
            cellID = args.input10X
            checkInputs(inputFolder, outputFolder, sampleSheet, cellID)
            processFastqFiles(inputFolder, outputFolder, sampleSheet, outputFileName, cellID)
    else:
        parser.print_help()
