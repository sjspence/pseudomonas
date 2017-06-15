#!/usr/bin/env python

import os
from optparse import OptionParser
from Bio import SeqIO
from Bio import Entrez
import numpy as np
import pandas as pd

def dirSlash(directory):
    if directory[-1] != '/':
	directory = directory + '/'
    return directory

def getFeatureName(seqFeature):
    quals = seqFeature.qualifiers
    name = np.nan
    if 'gene' in quals:
        name = quals['gene'][0]
    elif 'product' in quals:
        name = quals['product'][0]
    elif 'locus_tag' in quals:
        name = quals['locus_tag'][0]
    else:
        print('ERROR: missed feature')
        print(quals)
    return name

def recoverFeatures(gbFile, pos):
    records = list(SeqIO.parse(gbFile, 'genbank'))
    pos = int(pos)
    passedMapping = False
    for r in records:
        prevFeature = np.nan
        for f in r.features:
            if f.type == 'source':
                continue
            start = str(f.location).split(':')[0].replace('[', '')
            end = str(f.location).split(':')[1].split(']')[0]
            start = start.replace('>','').replace('<', '')   #Added replace since some locations are '>##'
            end = end.replace('>','').replace('<', '')
            start = start.replace('join{', '')
            end = end.replace('join{', '')
            strand = str(f.location).split('(')[1].replace(')', '')
            #Look for a coding feature
            if (pos > int(start)) and (pos < int(end)):
                name = getFeatureName(f)
                return {'coding': name}
            #Look for intergenic
            if (pos < int(start)) and (pos < int(end)):
                name = getFeatureName(f)
                prevName = np.nan
                if isinstance(prevFeature, float):
                    prevStrand = '+'
                else:
                    prevStrand = str(prevFeature.location).split('(')[1].replace(')', '')
                    prevName = getFeatureName(prevFeature)
                if (prevStrand == '+') and (strand == '+'):
                    return {'downstream1': prevName, 'upstream1': name}
                elif (prevStrand == '-') and (strand == '-'):
                    return {'downstream1': name, 'upstream1': prevName}
                elif (prevStrand == '+') and (strand == '-'):
                    return {'upstream1': prevName, 'upstream2': name}
                else:
                    return {'downstream1': prevName, 'downstream2': name}
            #store prev
            prevFeature = f

####
#Read in BLAST file
def blastToFeature(blastFile, gbDir, outFile):
    featureDf = pd.DataFrame()
    featureDict = {}
    f = open(blastFile, 'r')
    i = 0
    j = 1
    for line in f:
        if (i % 1000 == 0):
            print('Finished ' + str(j*i) + ' records.')
            j += 1
        i += 1
        line = line.strip().split('\t')
        readID = line[0]
        if readID in featureDict:
            continue
        gb = line[1].split('|')[3].split('.')[0]
        refStart = line[8]
        refEnd = line[9]
        gbFile = gbDir + gb + '.gb'
        featureDict[readID] = recoverFeatures(gbFile, refStart)
    df = pd.DataFrame.from_dict(featureDict)
    df.to_csv(outFile, sep = '\t')

def main():
    parser = OptionParser()
    parser.add_option("-i", "--inputFile", dest="i",
                      help="Input BLAST tabular file.",
                      metavar="INFILE")
    parser.add_option("-g", "--genbankDir", dest="g",
                      help="Directory with reference genbank files.",
                      metavar="GENBANK")
    parser.add_option("-o", "--outputFile", dest="o",
		      help="Output file to summarize gene mappings.",
		      metavar="OUTFILE")
    (options, args) = parser.parse_args()
    blastToFeature(options.i, options.g, options.o)

if __name__ == '__main__':
    main()
