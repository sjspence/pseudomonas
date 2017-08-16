#!/usr/bin/env python

from optparse import OptionParser
import pandas as pd

def downsampleSam(samFileName, diffFileName, outFileName):
    #Prep differentially mapped reads to look for
    readFwdRevSet = set()
    diffFile = open(diffFileName, 'r')
    for line in diffFile:
	readFwdRevSet.add(line.strip())
    diffFile.close()
    readSet = set()
    for i in readFwdRevSet:
        readSet.add(i.split('__')[0])
    #Search for those reads in a sam file
    samFile = open(samFileName, 'r')
    outFile = open(outFileName, 'w')
    for line in samFile:
        if '@SQ' == line[0:3]:
            outFile.write(line)
            continue
        lineList = line.strip().split('\t')
        if lineList[0] not in readSet:
            continue
        else:
            read = lineList[0]
            flag = int(lineList[1])
            if (flag & 64) and (read + '__1' in readFwdRevSet): #fwd read
                outFile.write(line)
            if (flag & 128) and (read + '__2' in readFwdRevSet): #rev read
                outFile.write(line)
    samFile.close()
    outFile.close()

def main():
    parser = OptionParser()
    parser.add_option("-s", "--samfile", dest="s",
                      help="SAM file.",
                      metavar="SAM")
    parser.add_option("-d", "--difffile", dest="d",
		      help="File with differentially mapped reads.",
		      metavar="DIFF")
    parser.add_option("-o", "--output", dest="o",
                      help="Output filtered SAM file.",
                      metavar="OUTPUT")
    (options, args) = parser.parse_args()
    downsampleSam(options.s, options.d, options.o)

if __name__ == '__main__':
    main()
