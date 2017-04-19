#!/usr/bin/env python

import os
from Bio import SeqIO
from Bio import AlignIO

def runAmphora2(inDir, checkDir, outDir, parallel):
    cleanDir = '08_all_clean_contigs/'
    checkFile = open(checkDir + 'summary.txt', 'r')
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
    if not os.path.exists(cleanDir):
        os.makedirs(cleanDir)
        for f in os.listdir(inDir):
            sampID = f.replace('_contigs.fna', '')
            if sampID not in badGenomes:
                os.system('cp ' + inDir + f + ' ' + cleanDir + f)
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outFileName = outDir.replace('/', '.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Amphora2 commands\n\n')
    #outFile.write('trap \"kill 0\" SIGINT\n')
    cwd = os.getcwd() + '/'
    for f in os.listdir(cleanDir):
        sampleOut = outDir + f.replace('_contigs.fna', '/')
        if not os.path.exists(sampleOut):
            os.makedirs(sampleOut)
        else:
            os.system('rm -r ' + sampleOut)
            os.makedirs(sampleOut)
        try:
            os.system('rm ' + cwd + cleanDir + f + '.orf')
        except OSError:
            pass
        outFile.write('sh -c \'cd ' + sampleOut + ' && exec ')
        outFile.write('MarkerScanner.pl -Bacteria -DNA ' + cwd + cleanDir + \
                        f + '\'\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + str(parallel))

def runAmphora2Align(inDir, parallel):
    outFileName = inDir.replace('/', '_align.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Amphora2 MarkerAlignTrim commands\n\n')
    cwd = os.getcwd() + '/'
    for subdir, dirs, files in os.walk(inDir):
        for d in dirs:
            cdDir = inDir + d + '/'
            outFile.write('sh -c \'cd ' + cdDir + ' && exec ')
            outFile.write('MarkerAlignTrim.pl -OutputFormat fasta\'\n')
	    #outFile.write('MarkerAlignTrim.pl\'\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + str(parallel))

def getAmphoraDNA(contigDir, amphoraDir, parallel):
    outFileName = amphoraDir.replace('/', '_DNA.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Amphora2 DNA sequence recovery commands\n\n')
    for subdir, dirs, files in os.walk(amphoraDir):
        for d in dirs:
	    sampleID = d
	    contigFile = contigDir + d + '_contigs.fna'
	    outFile.write('./scripts/amphora_getDNA.py -a ' + \
			  amphoraDir + d + ' -c ' + contigFile + '\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel --verbose -j ' + \
	      str(parallel))

def alignAmphoraDNA(amphoraDir, parallel):
    outFileName = amphoraDir.replace('/', '_DNAalign.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Amphora2 DNA sequence alignment commands\n\n')
    fxnalGenes = {}
    allCt = 0
    for subdir, dirs, files in os.walk(amphoraDir):
	for d in dirs:
	    allCt += 1
	    for f in os.listdir(amphoraDir + d):
		if '.fa' in f:
		    gene = f.replace('.fa', '')
		    if gene not in fxnalGenes:
			fxnalGenes[gene] = [d]
		    else:
			fxnalGenes[gene].append(d)
    for gene in fxnalGenes:
	#Only keep genes that exist in all samples
	if len(fxnalGenes[gene]) != allCt:
	    continue
	geneFile = open(amphoraDir + gene + '_all.fa', 'w')
	for d in fxnalGenes[gene]:
	    for line in open(amphoraDir + d + '/' + gene + '.fa', 'r'):
		if '>' in line:
		    geneFile.write(line.replace('>', '>' + d + ' '))
		else:
		    geneFile.write(line)
	geneFile.close()
    outFile.close()

def alignProtein(amphoraDir, parallel):
    outFileName = amphoraDir.replace('/', '_alignProt.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Amphora2 protein alignment commands\n\n')
    #Create multi-fasta for each marker gene that is present in all clean
    #genomes
    fxnalGenes = {}
    allCt = 0
    for subdir, dirs, files in os.walk(amphoraDir):
	for d in dirs:
	    allCt += 1
	    for f in os.listdir(amphoraDir + d):
		if '.pep' in f:
		    gene = f.replace('.pep', '')
		    if gene not in fxnalGenes:
			fxnalGenes[gene] = [d]
		    else:
			fxnalGenes[gene].append(d)
    for gene in fxnalGenes:
	if len(fxnalGenes[gene]) != allCt:
	    continue
	geneFile = open(amphoraDir + gene + '_pep.fa', 'w')
	#Consider that sometimes the files can have more than one seq
	records = list(SeqIO.parse(amphoraDir + gene + '_pep.fa', 'fasta'))
	if len(records) == 1:
	    SeqIO.write(records, amphoraDir + d + '/' + gene + '.pep', 'fasta')
	else:
	    for r in records:
		#START HERE, delete some below...
	for d in fxnalGenes[gene]:
	    for line in open(amphoraDir + d + '/' + gene + '.pep', 'r'):
		if '>' in line:
		    geneFile.write(line.replace('>', '>' + d + ' '))
		else:
		    geneFile.write(line)
	geneFile.close()
    #Align protein sequences with Muscle
    for f in os.listdir(amphoraDir):
	if '_pep.fa' in f:
	    outFile.write('muscle -in ' + amphoraDir + f + ' -fastaout ' + \
			  amphoraDir + f.replace('_pep.fa', '_pep_align.fa') + \
			  '\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel --verbose -j ' + \
              str(parallel))

def raxml(amphoraDir, parallel):
    outFileName = amphoraDir.replace('/', '_raxml.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Amphora2 raxml alignment commands\n\n')
    #Concatenate alignments
    dnaG = []
    for subdir, dirs, files in os.walk(amphoraDir):
	for d in dirs:
	    for f in os.listdir(amphoraDir + d):
		if 'dnaG.aln' in f:
		    alignments = AlignIO.read(open(amphoraDir + d + '/' + f, 'r'),
				'phylip-relaxed')
		    print(alignments.get_alignment_length())
		    dnaG.append(alignments)
    print(len(dnaG))
    SeqIO.write(alignments, 'trial.phy', 'phylip-relaxed')

def main():
    ###AMPHORA2
    inDir = '03_spades_q20_contigs/'
    checkDir = '04_checkm_spades_q20/'
    outDir = '09_clean_amphora2/'
    parallel = 38
    #runAmphora2(inDir, checkDir, outDir, parallel)
    #
    ###AMPHORA_ALIGN
    inDir = '09_clean_amphora2/'
    parallel = 38
    #runAmphora2Align(inDir, parallel)
    #
    ###GET_DNA
    contigDir = '08_all_clean_contigs/'
    amphoraDir = '09_clean_amphora2/'
    parallel = 38
    #getAmphoraDNA(contigDir, amphoraDir, parallel)
    #
    ###ALIGN_DNA
    amphoraDir = '09_clean_amphora2/'
    parallel = 38
    #alignAmphoraDNA(amphoraDir, parallel)
    #
    ###ALIGN_PROTEIN
    alignProtein(amphoraDir, parallel)
    #
    ###BUILD_TREE
    amphoraDir = '09_clean_amphora2/'
    parallel = 38
    #raxml(amphoraDir, parallel)

if __name__ == "__main__":
    main()
