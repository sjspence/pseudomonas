#!/usr/bin/env python

import os
import subprocess

def getSubgroups():
    subDict = {}
    subDict['subA'] = ['D17-102051', 'D17-102042', 'D17-102049', 'D17-102065','D17-102046',
	    'D17-102050', 'D17-102044', 'D17-102047', 'D17-102043','D17-102038',
	    'D17-102041', 'D17-102048', 'D17-102039', 'D17-102045','D17-102037',
	    'D17-102036', 'D17-102040']
    subDict['subB'] = ['D17-102052', 'D17-102058', 'D17-102074', 'D17-102197','D17-102209',
	    'D17-102095', 'D17-102141', 'D17-102054', 'D17-102195','D17-102194',
	    'D17-102057', 'D17-102205', 'D17-102085', 'D17-102061','D17-102094',
	    'D17-102202', 'D17-102208', 'D17-102237', 'D17-102207','D17-102206',
	    'D17-102211', 'D17-102210', 'D17-102072', 'D17-102092','D17-102086',
	    'D17-102158']
    subDict['subC'] = ['D17-102017', 'D17-102002', 'D17-102235','D17-102003', 'D17-102217']
    subDict['subD'] = ['D17-102232', 'D17-102007', 'D17-102228', 'D17-102008','D17-102006',
	    'D17-102176']
    subDict['subE'] = ['D17-102084', 'D17-102082', 'D17-102198','D17-102073', 'D17-102224']
    subDict['subF'] = ['D17-102076', 'D17-102033', 'D17-102060', 'D17-102014','D17-102075',
    	    'D17-102032', 'D17-102016', 'D17-102034', 'D17-102013','D17-102091',
    	    'D17-102035', 'D17-102214', 'D17-102015']
    subDict['subG'] = ['D17-102080', 'D17-102020', 'D17-102025', 'D17-102019','D17-102027',
	    'D17-102024', 'D17-102022', 'D17-102021', 'D17-102028','D17-102029',
	    'D17-102026', 'D17-102023']
    subDict['subH'] = ['D17-102088', 'D17-102090', 'D17-102215', 'D17-102089']
    subDict['subI'] = ['D17-102226', 'D17-102018', 'D17-102218', 'D17-102231','D17-102222']
    subDict['subJ'] = ['D17-102238', 'D17-102242', 'D17-102239', 'D17-102244','D17-102245',
	    'D17-102243', 'D17-102240', 'D17-102241', 'D17-102236','D17-102234']
    subDict['subK'] = ['D17-102004', 'D17-102180', 'D17-102220', 'D17-102069','D17-102070',
	    'D17-102212']
    subDict['subL'] = ['D17-102078', 'D17-102083', 'D17-102081', 'D17-102079']
    subDict['subM'] = ['D17-102203', 'D17-102055', 'D17-102204']
    subDict['subN'] = ['D17-102213', 'D17-102193', 'D17-102093']
    subDict['subO'] = ['D17-102010', 'D17-102009', 'D17-102229', 'D17-102230','D17-102005',
	    'D17-102233']
    return subDict

def setupDirs(s, samps):
    gutDir = '/scratch/users/mit_alm/gutevo/'
    sourceDir = gutDir + 'case_2017_05_11_D44M_SLN/'
    caseDir = gutDir + 'case_2017_06_16_sjs_' + s + '/'
    #1. Create case folder
    if not os.path.exists(caseDir):
	os.makedirs(caseDir)
    #2. Copy files
    os.system('cp ' + sourceDir + 'case.slurm ' + caseDir)
    #3. Change the genome name in build_mutation_table_master.m
    buildIn = sourceDir + 'build_mutation_table_master.m'
    inF = open(buildIn, 'r')
    buildOut = caseDir + 'build_mutation_table_master.m'
    outF = open(buildOut, 'w')
    genomeName = 'Pseudomonas_' + s
    refLine = 'REF_GENOME_DIRECTORY = [masterdir \'/Reference_Genomes/'
    for line in inF:
	if refLine in line:
	    line = refLine + genomeName + '/\'];\n'
	outF.write(line)
    inF.close()
    outF.close()
    #4. Change the sample names in sample_names.csv
    sampFile = caseDir + 'sample_names.csv'
    expDir = gutDir + '2017_06_16_sjs_' + s
    sampF = open(sampFile, 'w')
    sampF.write('ExperimentFolder,Sample,AlignmentFolder,Outgroup\n')
    for s in samps:
	sList = [expDir, s, s + '/sickle2050/' + genomeName, '0']
	sampF.write(','.join(sList) + '\n')
    sampF.close()
    #5. Submit job
    os.system('cd ' + caseDir + '\necho $PWD\nsbatch case.slurm')

def main():
    subDict = getSubgroups()
    for s in subDict:
	if (s == 'subA') or (s == 'subB') or (s == 'subG') or (s == 'subJ'):
	    continue
	setupDirs(s, subDict[s])

if __name__ == "__main__":
    main()
