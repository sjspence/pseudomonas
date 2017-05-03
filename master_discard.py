#!/usr/bin/env python

import os
import shutil
from Bio import SeqIO

def checkSlash(directory):
    if directory[-1] != '/':
        directory = directory + '/'
    return directory

def runAmphora2Align(inDir, parallel):
    outFileName = inDir.replace('/', '_align.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Amphora2 MarkerAlignTrim commands\n\n')
    cwd = os.getcwd() + '/'
    for subdir, dirs, files in os.walk(inDir):
        for d in dirs:
            cdDir = inDir + d + '/'
            outFile.write('sh -c \'cd ' + cdDir + ' && exec ')
            #outFile.write('MarkerAlignTrim.pl -OutputFormat fasta\'\n')
            outFile.write('MarkerAlignTrim.pl -WithReference -OutputFormat fasta\'\n')
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
#Edit the code below to consider that some genes grabbed smaller fragments too
#    fxnalGenes = {}
#    allCt = 0
#    for subdir, dirs, files in os.walk(amphoraDir):
#       for d in dirs:
#           allCt += 1
#           for f in os.listdir(amphoraDir + d):
#               if '.fa' in f:
#                   gene = f.replace('.fa', '')
#                   if gene not in fxnalGenes:
#                       fxnalGenes[gene] = [d]
#                   else:
#                       fxnalGenes[gene].append(d)
#    for gene in fxnalGenes:
#       #Only keep genes that exist in all samples
#       if len(fxnalGenes[gene]) != allCt:
#           continue
#       geneFile = open(amphoraDir + gene + '_dna.fa', 'w')
#       for d in fxnalGenes[gene]:
#           for line in open(amphoraDir + d + '/' + gene + '.fa', 'r'):
#               if '>' in line:
#                   geneFile.write(line.replace('>', '>' + d + ' '))
#               else:
#                   geneFile.write(line)
#       geneFile.close()

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

def main():
    ###AMPHORA_ALIGN
    inDir = '09_amphora2/'
    parallel = 38
    #runAmphora2Align(inDir, parallel)
    #
    ###GET_DNA
    contigDir = '08_all_clean_contigs/'
    amphoraDir = '09_amphora2/'
    parallel = 38
    getAmphoraDNA(contigDir, amphoraDir, parallel)
    #
    ###MASK_ALIGNMENTS
    inDir = '10_amphora_aligned/'
    #maskAlignments(inDir, parallel)
    #

if __name__ == "__main__":
    main()
