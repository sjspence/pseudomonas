#!/usr/bin/env python

import os
import shutil
from Bio import SeqIO
from Bio import AlignIO

def checkSlash(directory):
    if directory[-1] != '/':
	directory = directory + '/'
    return directory

#Runs the MarkerScanner Amphora2 script to recover marker genes
def runAmphora2(inDir, outDir, parallel):
    inDir = checkSlash(inDir)
    outDir = checkSlash(outDir)
    #Clean up from any past runs
    if os.path.exists(outDir):
	shutil.rmtree(outDir)
    for f in os.listdir(inDir):
	if f[-4:] == '.orf':
	    os.remove(inDir + f)
    #Create new directory for amphora results
    outFileName = outDir.replace('/', '.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Amphora2 commands\n\n')
    cwd = os.getcwd() + '/'
    os.makedirs(outDir)
    for f in os.listdir(inDir):
        sampleOut = outDir + f.replace('_contigs.fna', '/')
        if not os.path.exists(sampleOut):
            os.makedirs(sampleOut)
        outFile.write('sh -c \'cd ' + sampleOut + ' && exec ')
        outFile.write('MarkerScanner.pl -Bacteria -DNA ' + cwd + inDir + \
                        f + '\'\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + str(parallel))

def alignProtein(amphoraDir, outDir, parallel):
    #Create multi-fasta peptide file for each marker gene that is present in
    #all clean genomes
    fxnalGenes = {}
    allCt = 0
    amphoraDir = checkSlash(amphoraDir)
    #Collect all peptide files associated with an amphora gene
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
    #Check that the functional gene is in all genomes
    #Recover the longest representative sequence recovered by amphora2 for
    # each genome
    outDir = checkSlash(outDir)
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    for gene in fxnalGenes:
	if len(fxnalGenes[gene]) != allCt:
	    continue
	peptideList = []
	for d in fxnalGenes[gene]:
	    inPep = amphoraDir + d + '/' + gene + '.pep'
	    records = list(SeqIO.parse(inPep, 'fasta'))
	    maxLen = 0
	    maxID = ''
	    for r in records:
		if len(r.seq) > maxLen:
		    maxLen = len(r.seq)
		    maxID = r.id
	    for r in records:
		if r.id == maxID:
		    r.id = d
		    r.description = d + ' ' + r.description
		    peptideList.append(r)
	outPep = outDir + gene + '_pep.fa'
	SeqIO.write(peptideList, outPep, 'fasta')
    #Align protein sequences with Muscle
    #Create executable file to run in parallel
    outFileName = amphoraDir.replace('/', '_alignProt.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Amphora2 protein alignment commands\n\n')
    for f in os.listdir(outDir):
	if '_pep.fa' in f:
	    outFile.write('muscle -in ' + outDir + f + ' -fastaout ' + \
			  outDir + f.replace('_pep.fa', '_pep_align.fa') + \
			  '\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel --verbose -j ' + \
              str(parallel))

#Mask alignments
def maskAlignments(inDir, parallel):
    inDir = checkSlash(inDir)
    outFileName = inDir.replace('/', '_mask.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Zorro masking commands\n\n')
    for f in os.listdir(inDir):
	if '_pep_align.fa' in f:
	    maskRaw = inDir + f.replace('_align.fa', '.mask_raw')
	    outFile.write('zorro ' + inDir + f  + ' > ' + maskRaw + '\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel --verbose -j ' + \
              str(parallel))

def formatAlignments(inDir, parallel):
    inDir = checkSlash(inDir)
    for f in os.listdir(inDir):
	pass

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
    ###CLEAN_CONTIGS
    inDir = '03_spades_q20_contigs/'
    checkMfile = '04_checkm_spades_q20/summary.txt'
    outDir = '08_clean_contigs/'
    #cleanContigs(inDir, checkMfile, outDir)
    #
    ###AMPHORA2
    inDir = '08_clean_contigs/'
    outDir = '09_amphora2/'
    parallel = 38
    #runAmphora2(inDir, outDir, parallel)
    #
    ###ALIGN_PROTEIN
    inDir = '09_amphora2/'
    outDir = '10_amphora_aligned/'
    #alignProtein(inDir, outDir, parallel)
    #
    ###MASK_ALIGNMENTS
    inDir = '10_amphora_aligned/'
    #maskAlignments(inDir, parallel)
    #
    ###FORMAT_ALIGNMENTS
    inDir = '10_amphora_aligned/'
    formatAlignments(inDir, parallel)
    #
    ###BUILD_TREE
    inDir = '09_amphora2/'
    parallel = 38
    #raxml(amphoraDir, parallel)

if __name__ == "__main__":
    main()
