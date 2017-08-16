#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
from optparse import OptionParser

#Returns pandas dataframe of variable read mappings for a single sample of source reads
#Input = mappedID: source reads for mapping are from this sample
#        unmapped: directory of unmapped results
#Columns = samples in subgroup, plus 'sequence_1' and 'sequence_2' for the fwd & rev sequences of the read
#Rows = Read IDs
#Data = one if even one read mapped, zero if both paired ends are unmapped
def binaryMapping(mappedID, unmapped):
    unmappedDict = {}
    unmappedSeqs = {}
    for f in os.listdir(unmapped):
        if 'mapped_' + mappedID in f:
            contig = f.split('_')[1]
            unmappedDict[contig] = {}
            inFile = open(unmapped + f, 'r')
            for line in inFile:
                line = line.strip().split('\t')
                readID = line[0]
                flag = int(line[1])
                seq = line[9]
                if (flag & 64) and (flag & 4):
                    unmappedSeqs[readID + '__1'] = seq
                    unmappedDict[contig][readID + '__1'] = 10
                elif (flag & 128) and (flag & 4):
                    unmappedSeqs[readID + '__2'] = seq
                    unmappedDict[contig][readID + '__2'] = 10
            inFile.close()
    df = pd.DataFrame.from_dict(unmappedDict)
    df = df.fillna(value=1)             #1 if mapped
    df = df.replace(10, 0)              #0 if unmapped
    ####
    #Remove reads that are unmapped (e.g. zeros) in all samples
    df = df.loc[(df != 0).any(axis = 1)]
    #Remove reads that are mapped (e.g. ones) in all samples
    df = df.loc[(df.sum(axis=1) != len(df.columns))]
    df['sequence'] = ''
    for i in list(df.index):
        df.loc[i, 'sequence'] = unmappedSeqs[i]
    return df

def main():
    parser = OptionParser()
    parser.add_option("-u", "--unmappedDir", dest="u",
                      help="Directory containing unmapped reads.",
                      metavar="UNMAPPED")
    parser.add_option("-m", "--mappedID", dest="m",
                      help="Source sample for mapped reads.",
                      metavar="MAPPED")
    parser.add_option("-b", "--binaryFile", dest="b",
		      help="Binary file of read mappings.",
		      metavar="BINARY")
    (options, args) = parser.parse_args()
    binaryFile = options.b
    if not os.path.exists(binaryFile):
	print('Starting ' + options.m)
        df = binaryMapping(options.m, options.u)
        df.to_csv(binaryFile, sep = '\t')
	print('Finished ' + options.m)

if __name__ == '__main__':
    main()
