#!/usr/bin/env python

import os
import sys

def runCutadapt(dataDir, parallel):
    parallel = str(parallel)
    if not os.path.exists('01_cutadapt/'):
        os.makedirs('01_cutadapt/')
    if '/' not in dataDir:
	dataDir += '/'
    outFile = open('01_cutadapt.sh', 'w')
    outFile.write('#!/bin/bash\n#Cutadapt commands\n\n')
    for subdir, dirs, files in os.walk(dataDir):
	for d in dirs:
	    sampName = d.replace('-2350G','')
	    for f in os.listdir(dataDir + d):
		if ('_1_' in f) and ('.fastq' in f):
		    f1 = dataDir + d + '/' + f
	    	if ('_2_' in f) and ('.fastq' in f):
		    f2 = dataDir + d + '/' + f
	    outFile.write('cutadapt -a CTGTCTCTTAT -A CTGTCTCTTAT ' + \
			'-o 01_cutadapt/' + sampName + '_1_trim.fastq ' + \
			'-p 01_cutadapt/' + sampName + '_2_trim.fastq ' + \
			f1 + ' ' + f2 + '\n\n')
    outFile.close()
    os.system('chmod +x 01_cutadapt.sh')
    os.system('cat 01_cutadapt.sh | parallel -j ' + parallel)

def runTrimmomatic(trimQ, threads, parallel):
    trimQ = str(trimQ)
    threads = str(threads)
    parallel = str(parallel)
    outDir = '02_trimmomatic_q' + trimQ + '/'
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    outFileName = '02_trimmomatic_q' + trimQ + '.sh'
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Trimmomatic commands\n\n')
    for f in os.listdir('01_cutadapt/'):
	if ('_1_' in f):
	    f = f.split('_')
	    sampName = f[0]
	    f1 = '01_cutadapt/' + sampName + '_1_trim.fastq'
	    f2 = '01_cutadapt/' + sampName + '_2_trim.fastq'
	    out1 = outDir + sampName + '_1_trimmomatic.fastq'
	    solo1 = out1.replace('.fastq', '_solo.fastq')
	    out2 = outDir + sampName + '_2_trimmomatic.fastq'
	    solo2 = out2.replace('.fastq', '_solo.fastq')
	    outFile.write('trimmomatic PE -threads ' + threads + \
			' -phred33 ' + f1 + ' ' + f2 + ' ' + \
			out1 + ' ' + solo1 + ' ' + out2 + ' ' + solo2 + ' ' + \
			'LEADING:3 TRAILING:3 SLIDINGWINDOW:5:' + trimQ + \
			' MINLEN:50\n\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + parallel)

def runSpades(inDir, threads, parallel):
    outDir = inDir.replace('02', '03').replace('trimmomatic', 'spades')
    threads = str(threads)
    parallel = str(parallel)
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    outFileName = outDir.replace('/', '.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Spades commands\n\n')
    for f in os.listdir(inDir):
	if ('_1_trimmomatic.fastq' in f):
	    f = f.split('_')
	    sampName = f[0]
	    f1 = inDir + sampName + '_1_trimmomatic.fastq'
	    f2 = inDir + sampName + '_2_trimmomatic.fastq'
	    outFile.write('spades.py -k 21,33,55,77 --careful -t ' + threads + \
			' -1 ' + f1 + ' -2 ' + f2 + \
			' -o ' + outDir + sampName + '\n\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + parallel)

def runCheckM(inDir, threads):
    threads = str(threads)
    #Move contigs to separate folder for checkM
    contigDir = inDir.replace('/', '_contigs/')
    if not os.path.exists(contigDir):
	os.makedirs(contigDir)
    for subdir, dirs, files in os.walk(inDir):
	for d in dirs:
	    if 'D17-' in d:
		f1 = inDir + d + '/contigs.fasta'
		o1 = contigDir + d + '_contigs.fna'
		os.system('cp ' + f1 + ' ' + o1)
    #Run checkm
    outDir = inDir.replace('03_', '04_checkm_')
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    os.system('checkm lineage_wf -t ' + threads + ' ' + contigDir + ' ' + \
		outDir)

def runQuast(inDir, threads, parallel):
    threads = str(threads)
    parallel = str(parallel)
    outDir = inDir.replace('03_', '04_quast_')
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    outFileName = outDir.replace('/', '.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Quast commands\n\n')
    for subdir, dirs, files in os.walk(inDir):
	for d in dirs:
	    if 'D17-' in d:
		f1 = inDir + d + '/contigs.fasta'
		o1 = outDir + d
		outFile.write('quast.py -L -t ' + threads + ' -o ' + o1 + \
				' ' + f1 + '\n\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + parallel)

def runRNAmmer(inDir, parallel):
    parallel = str(parallel)
    outDir = inDir.replace('03_', '05_rnammer_')
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    outFileName = outDir.replace('/', '.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#RNAmmer commands\n\n')
    for subdir, dirs, files in os.walk(inDir):
	for d in dirs:
	    if 'D17-' in d:
		f1 = inDir + d + '/contigs.fasta'
		outFile.write('rnammer -S bac -m ssu' + \
				' -gff ' + outDir + d + '.gff' + \
				' -xml ' + outDir + d + '.xml' + \
				' -h ' + outDir + d + '_hmm_report.txt' + \
				' -f ' + outDir + d + '_16s.fa ' + \
				f1 + '\n\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + parallel)

def checkSlash(directory):
    if directory[-1] != '/':
        directory = directory + '/'
    return directory

#Go through a contig directory and copy high completeness, low contamination
#contigs to a fresh directory
def cleanContigs(inDir, checkMfile, outDir):
    inDir = checkSlash(inDir)
    outDir = checkSlash(outDir)
    checkFile = open(checkMfile, 'r')
    badGenomes = set()
    header = True
    for line in checkFile:
        if header:
            header = False
            continue
        line = line.strip().split('\t')
        sampID = line[0].replace('_contigs', '')
        completeness = float(line[11])
        contamination = float(line[12])
        if (completeness < 95.0) or (contamination > 10.0):
            badGenomes.add(sampID)
    checkFile.close()
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        for f in os.listdir(inDir):
            sampID = f.replace('_contigs.fna', '')
            if sampID not in badGenomes:
                os.system('cp ' + inDir + f + ' ' + outDir + f)

def main():
    #dataDir = sys.argv[1]
    #
    ###CUTADAPT
    parallel = 30
    #runCutadapt(dataDir, parallel)
    #
    ###TRIMMOMATIC
    trimQ = 20
    threads = 4
    parallel = 9
    #runTrimmomatic(trimQ, threads, parallel)
    #
    ###SPADES
    inDir = '02_trimmomatic_q' + str(trimQ) + '/'
    threads = 12
    parallel = 3
    #runSpades(inDir, threads, parallel)
    #
    ###CHECKM
    inDir = '03_spades_q' + str(trimQ) + '/'
    threads = 20
    #runCheckM(inDir, threads)
    #
    ###QUAST
    inDir = '03_spades_q' + str(trimQ) + '/'
    threads = 4
    parallel = 9
    #runQuast(inDir, threads, parallel)
    #
    ###RNAMMER
    inDir = '03_spades_q' + str(trimQ) + '/'
    parallel = 38
    #runRNAmmer(inDir, parallel)
    #
    ###CLEAN_CONTIGS
    inDir = '03_spades_q' + str(trimQ) + '_contigs/'
    checkMfile = '04_checkm_spades_q20/summary.txt'
    outDir = '08_clean_contigs/'
    #cleanContigs(inDir, checkMfile, outDir)

if __name__ == "__main__":
    main()
