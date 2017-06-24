#!/usr/bin/env python

import os
import sys

def checkSlash(directory):
    if directory[-1] != '/':
        directory = directory + '/'
    return directory

def concatenateSubgroups(inDir, subDir, outDir, parallel):
    inDir = checkSlash(inDir)
    subDir = checkSlash(subDir)
    outDir = checkSlash(outDir)
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    bashName = outDir.replace('/', '_concat.sh')
    bashFile = open(bashName, 'w')
    bashFile.write('#!/bin/bash\n#Concatenation commands\n\n')
    for subdirs, dirs, files in os.walk(subDir):
	for d in dirs:
	    print('Joining fastq reads for ' + d + '...')
	    subIDs = []
	    for f in os.listdir(subDir + d):
		if f[-3:] == '.fa':
		    subIDs.append(f.replace('_contigs.fa', ''))
	    outF = outDir + d + '_1_trimmomatic.fastq'
	    outR = outDir + d + '_2_trimmomatic.fastq'
	    outFfiles = []
	    outRfiles = []
	    for s in subIDs:
		outFfiles.append(inDir + s + '_1_trimmomatic.fastq')
		outRfiles.append(inDir + s + '_2_trimmomatic.fastq')
	    bashFile.write('cat ' + ' '.join(outFfiles) + ' > ' + outF + \
				'\n')
	    bashFile.write('cat ' + ' '.join(outRfiles) + ' > ' + outR + \
				'\n')
    bashFile.close() 
    os.system('chmod +x ' + bashName)
    os.system('cat ' + bashName + ' | parallel -j ' + str(parallel))

def runSpades(inDir, outDir, threads, parallel):
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

def runProkka(genomeDir, outDir, parallel):
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outFileName = '09_prokka.sh'
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Prokka commands\n\n')
    for f in os.listdir(genomeDir):
        if f[-4:] != '.fna':
            continue
        if f.split('_')[0] in os.listdir(outDir):
            continue
        genome = f.replace('_contigs.fna', '')
        outProkka = outDir + genome + '/'
        outFile.write('prokka --outDir ' + outProkka + ' --prefix ' + \
                        genome + ' ' + genomeDir + f + '\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + str(parallel))

def main():
    ###MAKE SUBGROUP DIRECTORIES
    inDir = '02_trimmomatic_q20/'
    subDir = '10_subgroups/'
    outDir = '13_subgroup/'
    parallel = 30
    #concatenateSubgroups(inDir, subDir, outDir, parallel)
    #
    ###SPADES
    inDir = '13_subgroup/'
    outDir = '14_subgroup_assemblies/'
    threads = 2
    parallel = 15
    #runSpades(inDir, outDir, threads, parallel)
    #
    ###
    inDir = '14_subgroup_contigs/'
    outDir = '14_subgroup_prokka/'
    parallel = 20
    runProkka(inDir, outDir, parallel)
    #
    #
    #
    #
    #
    #
    #
    ###CHECKM
    #inDir = '03_spades_q' + str(trimQ) + '/'
    threads = 20
    #runCheckM(inDir, threads)
    #
    ###QUAST
    #inDir = '03_spades_q' + str(trimQ) + '/'
    threads = 4
    parallel = 9
    #runQuast(inDir, threads, parallel)
    #
    ###RNAMMER
    #inDir = '03_spades_q' + str(trimQ) + '/'
    parallel = 38
    #runRNAmmer(inDir, parallel)
    #
    ###CLEAN_CONTIGS
    #inDir = '03_spades_q' + str(trimQ) + '_contigs/'
    checkMfile = '04_checkm_spades_q20/summary.txt'
    outDir = '08_clean_contigs/'
    #cleanContigs(inDir, checkMfile, outDir)

if __name__ == "__main__":
    main()
