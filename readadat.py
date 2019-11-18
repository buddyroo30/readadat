#!/usr/bin/env python3

import csv
import re
import json
import argparse
import os.path
from os import path
import pandas as pd
import numpy as np

#Overview of SomaLogic adat files:
#There are 2 kinds of rows in the ^TABLE_BEGIN section: (1) column annotation rows (i.e. annotations about the proteins
#being measured, such as their ID, sequence, gene, Uniprot ID, etc.) and (2) sample annotation rows, containing
#both annotations about the samples themselves such as sample ID, sample name, plate location, etc. as well as
#all the SomaLogic readouts for all the proteins for that sample.
#Format of column annotation row:
#<empty sample annotation values tab separated, since no sample being annotated on this row>\t<Column Annotation Name>\t<column annotation values tab separated>
#Format of sample annotation row:
#<sample annotation values tab separated, for the sample in this row>\t<empty Column Annotation Name, since this row is giving the SomaLogic readout values>\t<SomaLogic readout values for each protein, tab separated>
#All the column annotation names/keys are given in the ^COL_DATA section
#and all the sample annotation names/keys are given in the ^ROW_DATA section
#The ^HEADER section gives simple key/value annotations (tab separated each on its own line)
#that provide info about the whole experiment.
#
#The return value of this is a Python Dict that contains these keys and values
#(as extracted from the adat file):
#Checksum - the checksum value from the adat value
#Metadata - The key/value pairs from the ^HEADER section
#SampleAnnotNames - a Python list giving the names of the sample annotation attributes
#SequenceAnnotNames - a Python list giving the names of the sequence/gene/protein annotation attributes
#SequenceData - Annotation values for the sequences/genes/proteins as a Pandas Dataframe
#SampleAndIntensityData - Sample annotations and their SOMAlogic readout values as a Pandas Dataframe
#The output was built to be similar to the output of the readat R package's readAdat function
def readAdat (adat_file, keepOnlyPasses = True, keepOnlySamples = True):

    adatVals = {}
    curSection = None
    sample_type_field_i = None
    row_check_field_i = None
    sectionTitles = {'^HEADER','^COL_DATA','^ROW_DATA','^TABLE_BEGIN'}
    adatVals['Metadata'] = {}

    sampleAndIntensityData = []
    sequenceData = []
    sequenceDataIndex = []

    colPassFlag = None
    
    csv_reader = csv.reader(open(adat_file, mode='r'),delimiter='\t')
    line_count = 0
    for row in csv_reader:
        if curSection == '^TABLE_BEGIN':
            rowDataLen = len(adatVals['SampleAnnotNames'])
            curColName = row[rowDataLen]
            if (curColName): #Column annotation row
                curColAnnotSlice = row[rowDataLen+1:]
                sequenceData.append(curColAnnotSlice)
                sequenceDataIndex.append(curColName)
                if keepOnlyPasses and curColName == 'ColCheck':
                    colPassFlag = [val == 'PASS' for val in curColAnnotSlice]
            else: #Sample annotation row, with its annotations and its Somalogic readout values for the genes/cols
                curSampleAnnotSlice = row[:rowDataLen]
                sampleAnnotNames = adatVals['SampleAnnotNames']
                if (curSampleAnnotSlice != sampleAnnotNames):
                    if keepOnlySamples and sample_type_field_i is not None:
                        curSampleType = curSampleAnnotSlice[sample_type_field_i]
                        if curSampleType != 'Sample':
                            continue
                    if keepOnlyPasses and row_check_field_i is not None:
                        curRowCheck = curSampleAnnotSlice[row_check_field_i]
                        if (curRowCheck != 'PASS'):
                            continue
                    curSomaReadoutVals = row[rowDataLen+1:]
                    if keepOnlyPasses and colPassFlag:
                        curSomaReadoutVals = [i for (i, v) in zip(curSomaReadoutVals, colPassFlag) if v]
                    sampleAndIntensityData.append(curSampleAnnotSlice + curSomaReadoutVals)
        elif row[0] in sectionTitles:
            curSection = row[0]
        elif curSection is None and row[0] == '!Checksum': #checksum value for the entire adat file you can use to check that it is not corrupted
            adatVals['Checksum'] = row[1]
        elif curSection == '^HEADER':
            adatVals['Metadata'][row[0]] = row[1]
        elif curSection in {'^ROW_DATA','^COL_DATA'}:
            nameOrType = row[0]
            row.pop(0)
            if curSection == '^ROW_DATA' and nameOrType == '!Name':
                sample_type_field_i = row.index('SampleType')
                row_check_field_i = row.index('RowCheck')
                adatVals['SampleAnnotNames'] = row

        line_count += 1

    if keepOnlyPasses and colPassFlag:
        for i in range(0,len(sequenceData)):
            sequenceData[i] = [seqDataVal for (seqDataVal, passFlag) in zip(sequenceData[i], colPassFlag) if passFlag]

    sequenceData = pd.DataFrame(sequenceData,index=sequenceDataIndex)
    adatVals['SequenceAnnotNames'] = sequenceDataIndex

    sampleAnnotNames = adatVals['SampleAnnotNames']
    seqIds = ["SeqId." + SeqIdVal for SeqIdVal in sequenceData.loc['SeqId']]
    sequenceData.columns = seqIds #assign SeqIds as column labels for sequenceData
    df_cols = sampleAnnotNames + seqIds
    sampleAndIntensityData = pd.DataFrame(sampleAndIntensityData,columns=df_cols)
    #reset the intensity value columns to be type numpy.float64
    resetTypes = {}
    for ri in range(0,len(seqIds)):
        resetTypes[seqIds[ri]] = 'float64'
    sampleAndIntensityData = sampleAndIntensityData.astype(resetTypes)

    adatVals['SequenceData'] = sequenceData
    adatVals['SampleAndIntensityData'] = sampleAndIntensityData

    return adatVals

def empty(val):
    if val == None or re.match(r"^\s*$", val):
        return True
    else:
        return False

def main ():

    #simple example, just reads in adat file and prints it

    adat_file = ""
    keepOnlyPasses = "True"
    keepOnlySamples = "True"

    parser = argparse.ArgumentParser()
    parser.add_argument('--adat_file', help="the path to the ADAT file to be read and parsed", type=str)
    parser.add_argument('--keepOnlyPasses', help="Keep only passes in the ADAT file (True or False)? default True", default="True",type=str)
    parser.add_argument('--keepOnlySamples', help="Keep only samples in the ADAT file (True or False)? default True", default="True",type=str)

    try:
        args = parser.parse_args()
        adat_file = args.adat_file.strip()
        keepOnlyPasses = args.keepOnlyPasses.strip()
        keepOnlySamples = args.keepOnlySamples.strip()
    except argparse.ArgumentError as arg_error:
        print(arg_error)
        raise Exception("Error in command line arguments in readadat.py\n")

    assert (not empty(adat_file) and path.exists(adat_file)), "For argument adat_file, you must provide a non-empty path to a file that exists"
    assert (not empty(keepOnlyPasses) and (keepOnlyPasses == "True" or keepOnlyPasses == "False")), "keepOnlyPasses must be True or False"
    assert (not empty(keepOnlySamples) and (keepOnlySamples == "True" or keepOnlySamples == "False")), "keepOnlySamples must be True or False"

    if (keepOnlyPasses == "True"):
        keepOnlyPasses = True
    else:
        keepOnlyPasses = False

    if (keepOnlySamples == "True"):
        keepOnlySamples = True
    else:
        keepOnlySamples = False
        

    adat = readAdat(adat_file, keepOnlyPasses = keepOnlyPasses, keepOnlySamples = keepOnlySamples)
    print(adat)

if __name__ == "__main__":
    main();   

