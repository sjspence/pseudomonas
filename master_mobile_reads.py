#!/usr/bin/env python

import os

def runUnmapped(inDir, outDir, parallel):
    print('Finding shared reads within subgroup: ' + inDir.split('/')[1])
    if not os.path.exists(outDir):
	os.makedirs(outDir)
    outFileName = '_'.join(outDir.split('/'))[:-1] + '.sh'
    outFile = open(outFileName, 'w')
    outFile.write('#!/bin/bash\n#Unmapped commands\n\n')
    for f in os.listdir(inDir):
	outFile.write('samtools view -f 4 ' + inDir + f + ' > ' + outDir + \
			f.replace('.sam', '_unmapped.sam') + '\n')
    outFile.close()
    os.system('chmod +x ' + outFileName)
    os.system('cat ' + outFileName + ' | parallel -j ' + str(parallel))

def main():
    ####
    #ASSESS UNMAPPED
    inDir = '11_bwa/subA/'
    outDir = '12_bwa_unmapped/subA/'
    parallel = 17
    runUnmapped(inDir, outDir, parallel)

if __name__ == '__main__':
    main()
