#!/usr/bin/env python3
#
# Author: Manuel Franco Sevilla
# Resave PIDCalib ntuples in 6.24 so that it can be opened in ROOT 5 

from subprocess import Popen, PIPE, STDOUT
import sys
import os
import os.path as op
from glob import glob
import shutil as sh

from argparse import ArgumentParser


#################################
# Command line arguments parser #
#################################

def parseInput():
    parser = ArgumentParser(
        description='Resave PIDCalib ntuples in 6.24 so that it can be opened in ROOT 5.'
    )

    parser.add_argument('-i', '--inFolder', default='/home/public/pidcalib_ntuples/remote/',
                        help='Folder with input ntuples.')
    parser.add_argument('-t', '--tag', default='', help='Tag to select subfolders with desired input ntuples.')
    parser.add_argument('-o', '--outFolder', default='/home/public/pidcalib_ntuples/remote/resaved',
                        help='Folder to store output ntuples.')

    return parser.parse_args()

    
## Add colors for terminal output
def cTerm(msg, color):
    num = 30
    if color == 'red':     num = 91
    if color == 'green':   num = 92
    if color == 'yellow':  num = 93
    if color == 'blue':    num = 94
    if color == 'magenta': num = 95
    if color == 'cyan':    num = 96
    return f'\033[{num};1m{msg}\033[0m'

## Run shell command and print output only to terminal
def runCmd(cmd):
    print('\n'+cTerm(' '.join(cmd),'magenta')+'\n')
    with Popen(cmd, stdout=PIPE, stderr=STDOUT, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    return p.returncode

## Create folder, remove if it already exists
def mkFolder(outFolder):
    if op.isdir(outFolder):
        print(cTerm(f'{outFolder} already exist! Removing it...','magenta'))
        sh.rmtree(outFolder, ignore_errors=True)
    os.makedirs(outFolder)
    

        
########
# Main #
########

if __name__ == '__main__':
    args = parseInput()

    # Create output folder
    if not op.isdir(args.outFolder):
        os.makedirs(args.outFolder)

    for folder in glob(args.inFolder+'/*'):
        if args.tag not in folder: continue

        subfolder = args.outFolder + '/' + op.basename(folder)
        mkFolder(subfolder)
        for ntp in glob(folder + '/*'):
            if '.root' not in ntp: continue
            runCmd(['./bin/resave_pidcalib_ntuple.exe', '-i', ntp, '-o', subfolder])
