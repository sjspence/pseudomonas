#!/usr/bin/env python

import os

def buildIndices(scaffoldDir, parallel):
    outFileName = scaffoldDir[0:len(scaffoldDir)-1] + '.sh'
    outFileName = outFileName.replace('/', '_')
    outFile = open(outFileName, 'w')
    for f in os.listdir(scaffoldDir):
	scaffold = scaffoldDir + f
	outFile.write('bwa index ' + scaffoldDir + f + '\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + str(parallel))

def runAligner(readDir, scaffoldDir, outDir, parallel):
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    contigList = []
    for f in os.listdir(scaffoldDir):
	if '.fa' in f[-3:]:
	    contigList.append(f.replace('_contigs.fa', ''))
    outFileName = outDir[0:len(outDir)-1] + '.sh'
    outFileName = outFileName.replace('/', '_')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Alignment commands\n\n')
    for sampReads in contigList:
	fwdReads = readDir + sampReads + '_1_trimmomatic.fastq'
	revReads = readDir + sampReads + '_2_trimmomatic.fastq'
	for sampScaffold in contigList:
	    scaffold = scaffoldDir + sampScaffold + '_contigs.fa'
	    outFile.write('bwa mem ' + scaffold + ' ' + fwdReads + ' ' + \
			revReads + ' | gzip -3 > ' + outDir + \
			'contig_' + sampScaffold + '_mapped_' + sampReads + \
			'.sam.gz\n')
	outFile.write('\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + str(parallel))

def runUnzip(bwaDir, parallel):
    outFileName = bwaDir[0:len(bwaDir)-1] + '_unzip.sh'
    outFileName = outFileName.replace('/', '_')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Unzipping commands\n\n')
    for f in os.listdir(bwaDir):
	outFile.write('gunzip ' + bwaDir + f + '\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + str(parallel))
    os.system('rm ' + outFileName)

def main():
    #####
    #IDENTIFY SUBGROUPS
    subDir = '10_subgroups/'
    subGroups = []
    for subdirs, dirs, files in os.walk(subDir):
	for d in dirs:
	    subGroups.append(d)
    subGroups.sort()
    #Adjust for what's already been run
    subGroups = subGroups[3:]
    #####
    #BUILD INDICES
    parallel = 2
    for s in subGroups:
        #scaffoldDir = '10_subgroups/subA/'
        #scaffoldDir = '10_subgroups/subB/'
        #scaffoldDir = '10_subgroups/subC/'
	scaffoldDir = subDir + s + '/'
        #buildIndices(scaffoldDir, parallel)
    ######
    #ALIGN
    readDir = '02_trimmomatic_q20/'
    bwaOut = '11_bwa/'
    parallel = 20
    for s in subGroups:
        #scaffoldDir = '10_subgroups/subA/'
        #outDir = '11_bwa/subA/'
        #scaffoldDir = '10_subgroups/subB/'
        #outDir = '11_bwa/subB/'
        #scaffoldDir = '10_subgroups/subC/'
        #outDir = '11_bwa/subC/'
	scaffoldDir = subDir + s + '/'
	outDir = bwaOut + s + '/'
	#print('Aligning reads to contigs in subgroup: ' + s)
        #runAligner(readDir, scaffoldDir, outDir, parallel)
    #####
    #UNZIP
    for s in subGroups:
	bwaDir =  bwaOut + s + '/'
	runUnzip(bwaDir, parallel)

if __name__ == '__main__':
    main()
