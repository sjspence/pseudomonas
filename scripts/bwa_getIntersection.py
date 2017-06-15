#!/usr/bin/env python

import os
from optparse import OptionParser

def intersectMapped(inDir, mappedID, outFileName):
    print('Checking out mapped reads from ' + mappedID + '...')
    intersectionSet = set()
    ct = 0
    for f in os.listdir(inDir):
	flagSet = set()
        if ('mapped_' + mappedID) in f:
            contigID = f.split('_')[1]
            contigSet = set()
            inFile = open(inDir + f, 'r')
            for line in inFile:
                if '@SQ\t' in line:
                    continue
                line = line.strip().split('\t')
                readID = line[0]
                contigSet.add(readID)
		flagSet.add(line[1])
	    for i in flagSet:
		print(i)
	    print(inFile)
	    print(len(contigSet))
	    print(list(contigSet)[0:10])
            inFile.close()
	    if ct == 0:
		intersectionSet = contigSet.copy()
	    else:
                intersectionSet = intersectionSet.intersection(contigSet)
	    ct += 1
	    print(len(intersectionSet))
    outFile = open(outFileName, 'w')
    for read in intersectionSet:
	outFile.write(read + '\n')
    outFile.close()

def main():
    parser = OptionParser()
    parser.add_option("-i", "--inDir", dest="i",
                      help="Input directory of alignments.",
                      metavar="IN")
    parser.add_option("-m", "--mappedID", dest="m",
                      help="Sample ID used for read mapping.",
                      metavar="MAPID")
    parser.add_option("-o", "--outFile", dest="o",
		      help="Out file of intersection reads.",
		      metavar="OUTFILE")
    (options, args) = parser.parse_args()
    intersectMapped(options.i, options.m, options.o)

if __name__ == '__main__':
    main()
