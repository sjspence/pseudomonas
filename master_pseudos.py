#!/usr/bin/env python

import os

def checkSlash(directory):
    if directory[-1] != '/':
        directory = directory + '/'
    return directory

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

def makeCuratedDir(subName, subList):
    newDir = '10_subgroups/'
    if not os.path.exists(newDir + subName):
        os.makedirs(newDir + subName)
    for f in subList:
        os.system('cp 08_pseudo_contigs/' + f + '_contigs.fa ' + newDir + \
                  subName + f + '_contigs.fa')

def curateDirectories():
    newDir = '10_subgroups/'
    subA = ['D17-102051', 'D17-102042', 'D17-102049', 'D17-102065','D17-102046',
	    'D17-102050', 'D17-102044', 'D17-102047', 'D17-102043','D17-102038',
	    'D17-102041', 'D17-102048', 'D17-102039', 'D17-102045','D17-102037',
	    'D17-102036', 'D17-102040']
    subB = ['D17-102052', 'D17-102058', 'D17-102074', 'D17-102197','D17-102209',
	    'D17-102095', 'D17-102141', 'D17-102054', 'D17-102195','D17-102194',
	    'D17-102057', 'D17-102205', 'D17-102085', 'D17-102061','D17-102094',
	    'D17-102202', 'D17-102208', 'D17-102237', 'D17-102207','D17-102206',
	    'D17-102211', 'D17-102210', 'D17-102072', 'D17-102092','D17-102086',
	    'D17-102158']
    subC = ['D17-102017', 'D17-102002', 'D17-102235','D17-102003', 'D17-102217']
    subD = ['D17-102232', 'D17-102007', 'D17-102228', 'D17-102008','D17-102006',
	    'D17-102176']
    subE = ['D17-102084', 'D17-102082', 'D17-102198','D17-102073', 'D17-102224']
    subF = ['D17-102076', 'D17-102033', 'D17-102060', 'D17-102014','D17-102075',
    	    'D17-102032', 'D17-102016', 'D17-102034', 'D17-102013','D17-102091',
    	    'D17-102035', 'D17-102214', 'D17-102015']
    subG = ['D17-102080', 'D17-102020', 'D17-102025', 'D17-102019','D17-102027',
	    'D17-102024', 'D17-102022', 'D17-102021', 'D17-102028','D17-102029',
	    'D17-102026', 'D17-102023']
    subH = ['D17-102088', 'D17-102090', 'D17-102215', 'D17-102089']
    subI = ['D17-102226', 'D17-102018', 'D17-102218', 'D17-102231','D17-102222']
    subJ = ['D17-102238', 'D17-102242', 'D17-102239', 'D17-102244','D17-102245',
	    'D17-102243', 'D17-102240', 'D17-102241', 'D17-102236','D17-102234']
    subK = ['D17-102004', 'D17-102180', 'D17-102220', 'D17-102069','D17-102070',
	    'D17-102212']
    subL = ['D17-102078', 'D17-102083', 'D17-102081', 'D17-102079']
    subM = ['D17-102203', 'D17-102055', 'D17-102204']
    subN = ['D17-102213', 'D17-102193', 'D17-102093']
    subO = ['D17-102010', 'D17-102009', 'D17-102229', 'D17-102230','D17-102005',
	    'D17-102233']
    makeCuratedDir('subA/', subA)
    makeCuratedDir('subB/', subB)
    makeCuratedDir('subC/', subC)
    makeCuratedDir('subD/', subD)
    makeCuratedDir('subE/', subE)
    makeCuratedDir('subF/', subF)
    makeCuratedDir('subG/', subG)
    makeCuratedDir('subH/', subH)
    makeCuratedDir('subI/', subI)
    makeCuratedDir('subJ/', subJ)
    makeCuratedDir('subK/', subK)
    makeCuratedDir('subL/', subL)
    makeCuratedDir('subM/', subM)
    makeCuratedDir('subN/', subN)
    makeCuratedDir('subO/', subO)

def XparsnpSub(subGroup, threads):
    #mumiMax = str(mumiMax)
    outDir = '/home/ubuntu/proc/sjspence/170105_PSE/11_parsnp_subgroups/' + \
	 subGroup + '/'
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    curatedDir = '10_subgroups/' + subGroup + '/'
#    fileDir = '/home/ubuntu/proc/sjspence/170105_PSE/10_pseudo_edited/'
#    if not os.path.exists(fileDir):
#	os.makedirs(fileDir)
#    genomeDir = '08_pseudo_contigs/'
#    for f in os.listdir(genomeDir):
#	outFile = open(fileDir + f, 'w')
#	for line in open(genomeDir + f, 'r'):
#	    if '>' in line:
#		outFile.write(line.replace('\n', ' PROKKA_20170405\n'))
#	    else:
#		outFile.write(line)
    #NOTE: had to edit reference to include prokka version & accession
    #      identifiers
    #genbankFile = reference + '.gbk'
    #os.system('parsnp -g ' + genbankFile + ' -d ' + fileDir + ' -o ' + \
#		outDir + ' -U 0.05 -p ' + str(threads))
    #os.system('parsnp -r ' + genomeDir + reference + '_contigs.fa -d ' + \
#		genomeDir + ' -U ' + mumiMax + ' -p ' + str(threads) + \
#		' -o ' + outDir)
    os.system('parsnp -r ! -d ' + \
		curatedDir + ' -c -p ' + str(threads) + \
		' -o ' + outDir)

def parsnpSub(inDir, outDir, threads, parallel):
    inDir = checkSlash(inDir)
    outDir = checkSlash(outDir)
    outFileName = inDir.replace('/', '.sh')
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Parsnp commands\n\n')
    for subdir, dirs, files in os.walk(inDir):
	for d in dirs:
	    subOutDir = outDir + d
	    if not os.path.exists(subOutDir):
		os.makedirs(subOutDir)
	    outFile.write('parsnp -r ! -d ' + inDir + d + ' -c -p ' + \
			str(threads) + ' -o ' + subOutDir + '\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + str(parallel))

def parsnpDefault(reference, inDir, outDir, threads):
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    os.system('parsnp -r ' + reference + ' -d ' + inDir + ' -p ' + threads + \
	      ' -o ' + outDir)

def main():
    #getPseudos()
    #runMugsy()
    #runProkka()
    #curateDirectories()
    #
    #ALIGN_GENOMES
    inDir = '10_subgroups'
    outDir = '11_parsnp_subgroups'
    threads = 3
    parallel = 13
    parsnpSub(inDir, outDir, threads, parallel)
    #
    reference = '08_pseudo_contigs/D17-102043_contigs.fa'
    inDir = '08_pseudo_contigs/'
    outDir = '10_parsnp_' + reference + '_default/'
    #parsnpDefault(reference, inDir, outDir, threads)

if __name__ == "__main__":
    main()
