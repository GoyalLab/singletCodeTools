import os
import glob
import pandas as pd
from Bio import SeqIO
import Levenshtein
import gzip

def checkInputs(inputFolder, outputFolder, sampleSheet, cellIDPath=None):
    # Check if input folder exists
    if not os.path.isdir(inputFolder):
        raise ValueError("Input folder does not exist: {}".format(inputFolder))
    
    # Check if output folder exists
    if not os.path.isdir(outputFolder):
        raise ValueError("Output folder does not exist: {}".format(outputFolder))
    
    # Check if sample sheet is a .csv file and exists
    if not sampleSheet.endswith('.csv') or not os.path.isfile(sampleSheet):
        raise ValueError("Sample sheet must be an existing .csv file.")
    
    # Read sample sheet and check for sampleName and sampleNumber columns
    sample_data = pd.read_csv(sampleSheet, header=0)
    if 'sampleName' not in sample_data.columns or 'sampleNumber' not in sample_data.columns:
        raise ValueError("Sample sheet must contain 'sampleName' and 'sampleNumber' columns.")
    
    # Check for the existence of both R1 and R2 .fastq.gz files for each sample
    for index, row in sample_data.iterrows():
        r1_filename_pattern = os.path.join(inputFolder, f"{row['sampleName']}_S{row['sampleNumber']}_L001_R1_001.fastq.gz")
        r2_filename_pattern = os.path.join(inputFolder, f"{row['sampleName']}_S{row['sampleNumber']}_L001_R2_001.fastq.gz")

        if len(glob.glob(r1_filename_pattern)) == 0:
            raise ValueError(f"Read 1 fastq file not found for sample {row['sampleName']} with number {row['sampleNumber']} in the expected format.")
        if len(glob.glob(r2_filename_pattern)) == 0:
            raise ValueError(f"Read 2 fastq file not found for sample {row['sampleName']} with number {row['sampleNumber']} in the expected format.")
    
    # Optional: Check cellID file if provided
    if cellIDPath is not None:
        if not os.path.isfile(cellIDPath):
            raise ValueError("Specified cellID file does not exist: {}".format(cellIDPath))
        cellIDFileName = os.path.basename(cellIDPath)
        if cellIDFileName not in ['barcodes.tsv', 'barcodes.tsv.gz']:
            raise ValueError("cellID file must be named 'barcodes.tsv' or 'barcodes.tsv.gz'.")
    
    print("All the inputs for the command are valid and will proceed with creating the barcode sheet for all the samples in the sheet.")

def processFastqFiles(inputFolder, outputFolder, sampleSheetPath, outputFileName, cellID = None):
    sampleSheet = pd.read_csv(sampleSheetPath)
    if(cellID):
        validBarcodes10X = set(pd.read_csv(cellID, header=None)[0])
    else:
        validBarcodes10X = None
    allSamplesDataFrame = pd.DataFrame()
    # Process each sample
    for index, row in sampleSheet.iterrows():
        barcodeSample = processSampleBarcode(row['sampleName'], row['sampleNumber'], 
                                      inputFolder, validBarcodes10X)
        allSamplesDataFrame = pd.concat([allSamplesDataFrame, barcodeSample]
                                        , ignore_index=True)
    
    finalDF = allSamplesDataFrame.loc[allSamplesDataFrame.index.repeat(allSamplesDataFrame['count'])]
    finalDF = finalDF.drop(columns=['count'])
    if not outputFileName.endswith('.csv'):
        outputFileName += '.csv'
    
    finalDF.to_csv(os.path.join(outputFolder, outputFileName),
                   index=False)




def processSampleBarcode(sampleName, sampleNumber, inputFolder, cellIDList=None):
    r1FilePath = os.path.join(inputFolder, f"{sampleName}_S{sampleNumber}_L001_R1_001.fastq.gz")
    r2FilePath = os.path.join(inputFolder, f"{sampleName}_S{sampleNumber}_L001_R2_001.fastq.gz")
    
    # Initialize a list to hold parsed data
    read1Seq = []
    read2Seq = []

    # Open and process the gzip-compressed FASTQ files
    with gzip.open(r1FilePath, 'rt') as r1File:  # 'rt' mode for text reading
        for record in SeqIO.parse(r1File, "fastq"):
            read1Seq.append(str(record.seq))

    with gzip.open(r2FilePath, 'rt') as r2File:
        for record in SeqIO.parse(r2File, "fastq"):
            read2Seq.append(str(record.seq))
    df = pd.DataFrame({'Read1': read1Seq, 'Read2': read2Seq})
    #extracting cell ID from read 1
    df['cellID'] = df['Read1'].str.extract(r'^([ATCG]{16})')
    df['cellID'] = df['cellID'].astype(str) + "-" + str(sampleNumber)
    #keeping only unique rows and thus removing any read information
    dfUnique = df.drop_duplicates()
    #extracting UMI from read 1 and making sure it is valid
    if(cellIDList):
        dfFiltered = dfUnique[dfUnique['cellID'].isin(cellIDList)].copy()
    else: 
        dfFiltered = dfUnique.copy()
    print("Filtered rows of dataframe:", len(dfFiltered))
    dfFiltered['umi'] = dfFiltered['Read1'].str[-10:]
    dfFiltered = dfFiltered[dfFiltered['umi'].str.match(r'^[ACTG]{10}$')]
    #Identifying lineage barcode information and 
    pattern = r"GGGCTG(([AT][CG]|[CG][AT]){15})GACGCT"
    dfFiltered['barcode'] = dfFiltered['Read2'].str.extract(pattern)[0]
    dfFiltered['barcode'] = dfFiltered['barcode'].astype(str)

    # Filter rows based on the length of the barcode
    # Since we're looking for a specific length (30), we ensure that we only keep rows meeting this criteria.
    dfFiltered =dfFiltered[dfFiltered['barcode'].apply(lambda x: len(x) == 30)]

    dialOut = dfFiltered.groupby(['cellID', 'barcode']).agg(count =('umi', 'nunique')).reset_index()
    # dialOut = dialOut.groupby('cellID').apply(lambda x: collapseLineageBarcodes(x[['barcode', 'count']])).reset_index(drop=True)

    updatedColumns = dialOut.groupby('cellID').apply(lambda x: collapseLineageBarcodes(x[['barcode', 'count']])).reset_index(drop=True)

    # If the transformed_cols is not aligned (same index) with your original DataFrame, you might need to adjust it.
    # Assuming transformed_cols now contains transformed 'barcode' and 'count' with the same index as dialOut_2.

    # Update the original DataFrame with the transformed values
    dialOut.update(updatedColumns)

    dialOut["sample"] = sampleName
    finalData = dialOut[['cellID', 'barcode', 'count', 'sample']].copy()
    return finalData


def collapseLineageBarcodes(single10xBarcodeDf):
    if single10xBarcodeDf.empty:
        print("Input data frame is empty.")
        return single10xBarcodeDf
    
    # How many UMI support this 10x barcode
    originalUmiCount = single10xBarcodeDf['count'].sum()
    
    # Order such that the most supported pair would be the first
    single10xBarcodeDf = single10xBarcodeDf.sort_values(by='count', ascending=False)
    
    # Levenshtein Distance of all lineage barcodes
    distances = [Levenshtein.distance(bc, single10xBarcodeDf.iloc[0]['barcode']) for bc in single10xBarcodeDf['barcode']]
    single10xBarcodeDf['edistToMostCommonLB'] = distances
    
    # Remove any row that has an edit distance of 3 or less
    single10xBarcodeDf = single10xBarcodeDf[~((single10xBarcodeDf['edistToMostCommonLB'] < 4) & (single10xBarcodeDf['edistToMostCommonLB'] > 0))]
    
    # Calculate how many UMI were subtracted and add these UMIs to the most common barcode
    single10xBarcodeDf.loc[single10xBarcodeDf.index[0], 'count'] += originalUmiCount - single10xBarcodeDf['count'].sum()
    
    # If after collapsing the less supported lineages are supported only by one barcode and the leading barcode is larger than 5-remove these rows
    topCount = single10xBarcodeDf.iloc[0]['count']
    single10xBarcodeDf = single10xBarcodeDf[~((single10xBarcodeDf['count'] * 3) < topCount)]
    return single10xBarcodeDf



