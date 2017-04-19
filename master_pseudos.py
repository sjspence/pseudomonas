#!/usr/bin/env python

import os

def getPseudos():
    #check for really bad genomes
    badGenomes = set()
    checkFile = open('04_checkm_spades_q20/summary.txt', 'r')
    ct = 0
    for line in checkFile:
	if ct == 0:
	    ct += 1
	    continue
	line = line.strip().split('\t')
	sampID = line[0].replace('_contigs', '')
	completeness = float(line[11])
	contamination = float(line[12])
	if (completeness < 95.0) or (contamination > 10.0):
	    badGenomes.add(sampID)
	ct += 1
    checkFile.close()
    dir16s = '06_consensus_16s/'
    for f in os.listdir(dir16s):
	if '.sintax' in f:
	    with open(dir16s + f, 'r') as t:
		for line in t:
		    line = line.strip().split('\t')
		    sampID = f.replace('.sintax', '')
		    if ('g:Pseudomonas' in line[3]) and \
			(sampID not in badGenomes):
			os.system('cp 03_spades_q20/' + sampID + \
			'/contigs.fasta 08_pseudo_contigs/' + sampID + \
			'_contigs.fa')

def runMugsy():
    outDir = '/home/ubuntu/proc/sjspence/170105_PSE/09_mugsy/'
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    genomeDir = '08_pseudo_contigs/'
    genomes = []
    for f in os.listdir(genomeDir):
	genomes.append(genomeDir + f)
    outFile = open('09_mugsy.sh', 'w')
    outFile.write('#!/bin/bash\n#Mugsy commands\n\n')
    outFile.write('mugsy --directory ' + outDir + ' --prefix pseudo_genomes ' \
		+ ' '.join(genomes))
    outFile.close()
    os.system('chmod +x 09_mugsy.sh')
    #os.system('bash')
    #os.system('./home/ubuntu/tools/mugsy_x86-64-v1r2.3/mugsyenv.sh')
    os.system('./09_mugsy.sh')

def runProkka():
    genomeDir = '08_pseudo_contigs/'
    outDir = '09_prokka/'
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    outFileName = '09_prokka.sh'
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Prokka commands\n\n')
    for f in os.listdir(genomeDir):
	genome = f.replace('_contigs.fa','')
	outProkka = outDir + genome + '/'
	outFile.write('prokka --outDir ' + outProkka + ' --prefix ' + \
			genome + ' ' + genomeDir + f + '\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('./' + outFileName)

def runParsnp(reference, threads):
    mumiMax = str(0.50)
    outDir = '/home/ubuntu/proc/sjspence/170105_PSE/10_parsnp_' + reference + \
		'_' + mumiMax + '/'
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    fileDir = '/home/ubuntu/proc/sjspence/170105_PSE/10_pseudo_edited/'
    if not os.path.exists(fileDir):
	os.makedirs(fileDir)
    genomeDir = '08_pseudo_contigs/'
    for f in os.listdir(genomeDir):
	outFile = open(fileDir + f, 'w')
	for line in open(genomeDir + f, 'r'):
	    if '>' in line:
		outFile.write(line.replace('\n', ' PROKKA_20170405\n'))
	    else:
		outFile.write(line)
    #NOTE: had to edit reference to include prokka version & accession
    #      identifiers
    genbankFile = reference + '.gbk'
    #os.system('parsnp -g ' + genbankFile + ' -d ' + fileDir + ' -o ' + \
#		outDir + ' -U 0.05 -p ' + str(threads))
    os.system('parsnp -r ' + genomeDir + reference + '_contigs.fa -d ' + \
		genomeDir + ' -U ' + mumiMax + ' -p ' + str(threads) + \
		outDir)

def main():
    #getPseudos()
    #runMugsy()
    #runProkka()
    threads = 38
    runParsnp('D17-102043', 38)

if __name__ == "__main__":
    main()
