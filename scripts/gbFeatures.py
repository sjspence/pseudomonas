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
    name = ''
    nameList = []
    for q in quals:
        if q == 'translation':
            continue
        qTextList = quals[q]
        qTextString = '--'.join(qTextList)
        nameList.append(q + '||' + qTextString)
    name = '__'.join(nameList)
    return name

def recoverFeatures(readID, gbFile, qStart, qEnd, o):
    records = list(SeqIO.parse(gbFile, 'genbank'))
    qStart = int(qStart)
    qEnd = int(qEnd)
    coding = False
    types = set()
    for r in records:
        prevFeature = np.nan
        for i, f in enumerate(r.features):
            #Check for coding
            if (qStart in f) or (qEnd in f):
                if (f.type != 'source') and (f.type != 'CDS'):
                    types.add(f.type)
                    continue
                outList = [readID, qStart, qEnd, f.type, f.location.start, f.location.end, f.location.strand, 
                          getFeatureName(f)]
                outList = map(str, outList)
                o.write('\t'.join(outList) + '\n')
            #Check for up/downstream???
    return o, types

def annotateBlast(blastFile, gbDir, outFile):
    f = open(blastFile, 'r')
    o = open(outFile, 'w')
    reads = set()
    typesAll = set()
    for line in f:
        line = line.strip().split('\t')
        readID = line[0]
        if readID in reads: #Keep only top BLAST hit
            continue
        reads.add(readID)
        gb = line[1].split('|')[3].split('.')[0]
        refStart = line[8]
        refEnd = line[9]
        gbFile = gbDir + gb + '.gb'
        o, types = recoverFeatures(readID, gbFile, refStart, refEnd, o)
        for t in types:
            typesAll.add(t)
    print(typesAll)
    f.close()
    o.close()

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
    annotateBlast(options.i, options.g, options.o)

if __name__ == '__main__':
    main()
