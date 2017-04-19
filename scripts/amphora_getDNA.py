#!/usr/bin/env python

import os
from optparse import OptionParser

compDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def dirSlash(directory):
    if directory[-1] != '/':
	directory = directory + '/'
    return directory

def revComplement(seq):
    compSeq = ''
    for letter in seq:
	newLetter = compDict[letter]
	compSeq += newLetter
    return compSeq[::-1]

def importFasta(fileName):
    reads = {}
    with open(fileName, 'r') as f:
	header = ''
	seq = ''
	for line in f:
	    if '>' in line:
		if header != '':
		    reads[header] = seq
		header = line.strip()
		seq = ''
	    else:
		seq += line.strip()
    return reads

def getAmphoraDNA(amphoraDir, contigFile):
    amphoraDir = dirSlash(amphoraDir)
    for f in os.listdir(amphoraDir):
	if '.pep' in f:
    	    outFile = open((amphoraDir + f).replace('.pep', '.fa'), 'w')
    	    header, contig, start, end = '', '', '', ''
	    reverse = False
    	    with open(amphoraDir + f, 'r') as pepFile:
    	        for line in pepFile:
    	            if '>' not in line:
    	        	continue
    	            header = line
    	            line = line.strip().split(' ')
    	            contig = line[0]
		    #Reverse Sense
		    if (len(line) == 6) and ('REVERSE' in line[4]):
			reverse = True
    	            start = int(line[1].replace('[', ''))
    	            end = int(line[3].replace(']', ''))
	    finalContigs = importFasta(contigFile)
	    for h in finalContigs:
		if h in contig:
		    outFile.write(header)
		    if not reverse:
    	        	outFile.write(finalContigs[h][start-1:end] + '\n')
		    else:
			sourceDNA = revComplement(finalContigs[h][end-1:start])
			outFile.write(sourceDNA + '\n')
    	    outFile.close()

def main():
    parser = OptionParser()
    parser.add_option("-a", "--amphoraDir", dest="a",
                      help="Amphora peptide directory.",
                      metavar="AMPHORA")
    parser.add_option("-c", "--contigFile", dest="c",
                      help="Contigs to search for amphora seqs.",
                      metavar="CONTIGS")
    (options, args) = parser.parse_args()
    getAmphoraDNA(options.a, options.c)

if __name__ == '__main__':
    main()
