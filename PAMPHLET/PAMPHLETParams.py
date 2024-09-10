#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
PAMPHLET - PAM Prediction HomoLogous Enhancement Toolkit
'''

import argparse
import sys
import os
from Bio import SeqIO

def get_arg_parser():

    ap = argparse.ArgumentParser(description='PAMPHLET - PAM Prediction HomoLogous Enhancement Toolkit')

    '''
    ap.add_argument("-s","--spacer",help="Spacer sequences, required.",required=True)
    ap.add_argument("-r","--repeat",help="Repeat sequence, not required. Repeat sequence will using to revise flank sequence.",default=None)
    ap.add_argument("-p","--protein",help="Related class II Cas effector protein sequence, not required. Need to co-operate with --repeat. Need to set as a FASTA FILE.",default=None)
    ap.add_argument("-o","--outdir",help="Output directory, required.",required=True)
    ap.add_argument("-O","--orientation",help="Orientation of the repeat sequence, choice=negative/positive. Default is [positive].",choices=["negative","positive"],default="positive")
    ap.add_argument("-u","--unique",help="Unique mode. If set, only revise unique spacer. Default is [False].",action="store_true",default=False)
    ap.add_argument("-P","--proteinblast",help="Use this parameters if you already have outfmt 6 protein blast result(vs. taxid:2).")
    ap.add_argument("-l","--proteinlen",help="Length of protein sequence, required if --proteinblast is set.",type=int,default=None)
    ap.add_argument("-d","--spacerdb",help="Spacer blast database, default is [phage].",choices=['phage','prokaryote'],default='phage')
    ap.add_argument("-L","--flanklen",help="Length of flank sequence, default is [12].",type=int,default=12)
    ap.add_argument("-R","--refmode",help="Reference mode. If set, revise protein homolog based on taxid. r is reference and nr is non-reference. Default is [r]",choices=["r","nr","a"],default="r")
    ap.add_argument("-f","--freqmode",help="Base frequency calculation mode. Default is [sigmoid].",choices=['sigmoid','linear'],default='sigmoid')
    ap.add_argument("-b","--blastmode",help="Spacer blastn mode. Common mode means use default blastn parameters and strict mode means use specific parameters, which could get more putative protospacers but also could cause false positives. Default is [common].",choices=['relax','common'],default='common')
    ap.add_argument("--pcovs",help="Minimum percent coverage of spacer sequence, default is [0.9].",type=float,default=0.9)
    ap.add_argument("--pident",help="Minimum percent identity of spacer sequence, default is [0.9].",type=float,default=0.9)
    ap.add_argument("--rident",help="Minimum percent identity of repeat sequence, default is [0.8].",type=float,default=0.8)
    ap.add_argument("--MaxProteinNum",help="Maximum number of protein homologs, default is [20]. If the size is too large, NCBI will forbidden your request.",type=int,default=20)
    ap.add_argument("--logfile",help="Generate a logfile.",action="store_true")
    ap.add_argument("--reviseLen",help="Revise length of spacer sequence, default is [FALSE].",action="store_true")
    ap.add_argument("--quiet",help="Quiet mode. If set, do not print log to screen.",action="store_true")
    '''

    ap.add_argument("-s","--spacer",help="Spacer sequences, required.",required=True)
    ap.add_argument("-r","--repeat",help="Repeat sequence. Repeat sequence will using to revise flank sequence.",default=None, required=True)
    ap.add_argument("-p","--protein",help="Related class II Cas effector protein sequence, not required. Need to co-operate with --repeat. Need to set as a FASTA FILE.",default=None,required=True)
    ap.add_argument("-o","--outdir",help="Output directory, required.",required=True)
    ap.add_argument("-bDB",'--bacteriaDB',help="Bacteria database for PAMPHLET-LOCAL, for homologous CRISPR-Cas search and spacer expand.",required=True)
    ap.add_argument("-pDB",'--protoDB',help="Protospacer source genome database for PAMPHLET-LOCAL.",default="")
    ap.add_argument("-cDB",'--casDB',help="Cas protein database for PAMPHLET-LOCAL.",default="")
    ap.add_argument("-L","--flanklen",help="Length of flank sequence, default is [12].",type=int,default=12)
    ap.add_argument("-b","--blastmode",help="Spacer blastn mode. Common mode means use default blastn parameters and strict mode means use specific parameters, which could get more putative protospacers but also could cause false positives. Default is [common].",choices=['relax','common'],default='common')
    ap.add_argument("--pcovs",help="Minimum percent coverage of spacer sequence, default is [0.9].",type=float,default=0.9)
    ap.add_argument("--pident",help="Minimum percent identity of spacer sequence, default is [0.9].",type=float,default=0.9)
    ap.add_argument("--rident",help="Minimum percent identity of repeat sequence, default is [0.8].",type=float,default=0.8)
    ap.add_argument("--reviseLen",help="Revise length of spacer sequence, default is [FALSE].",action="store_true")
    ap.add_argument("-f","--freqmode",help="Base frequency calculation mode. Default is [sigmoid].",choices=['sigmoid','linear'],default='sigmoid')
    ap.add_argument("-O","--orientation",help="Orientation of the repeat sequence, choice=negative/positive. Default is [positive].",choices=["negative","positive"],default="positive")
    ap.add_argument("-bt","--BlastThreads",help="Number of threads for blastn, default is [1].",type=int,default=1)
    
    return ap

def check_spacer_sequences(spacer_file):
    illegal_bases = ['b','d','e','f','h','i','j','k','l','m','n','o','q','u','v','w','x','y','z','B','D','E','F','H','I','J','K','L','M','N','O','Q','U','V','W','X','Y','Z']
    for spacers in SeqIO.parse(spacer_file,'fasta'):
        for bases in str(spacers.seq):
            if bases in illegal_bases:
                print('Illegal base %s found in spacer sequence %s. Please check your spacer sequence file.' % (bases,spacers.id))
                sys.exit(1)

def check_environment():
    ### check minCED add to PATH ###
    if os.popen('which minced').read().strip() == '':
        print('minCED not found. Please add minced to PATH.')
        sys.exit(1)
    
    ### check seqkit add to PATH ###
    if os.popen('which seqkit').read().strip() == '':
        print('seqkit not found. Please add seqkit to PATH.')
        sys.exit(1)
    
    ### Check Weblogo add to PATH ###
    if os.popen('which weblogo').read().strip() == '':
        print('Weblogo not found. Please add weblogo to PATH.')
        sys.exit(1)

    if os.popen('which ghostscript').read().strip() == "" and os.popen('which gs').read().strip() == "":
        print("GhostScript not found. Please add GhostScript to PATH.")
        sys.exit(1)
        
    if os.popen('which blastp').read().strip() == '':
        print('blastp not found. Please add blastp to PATH.')
        sys.exit(1)
    
    if os.popen('which makeblastdb').read().strip() == '':
        print('makeblastdb not found. Please add makeblastdb to PATH.')
        sys.exit(1)
    
    if os.popen('which blastn').read().strip() == '':
        print('blastn not found. Please add blastn to PATH.')
        sys.exit(1)

def check_outdir(outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        print('Output directory %s already exists.' % outdir)
        sys.exit(1)

def check_protein_files(protein_file):
    seqNum = os.popen('grep -c ">" %s' % protein_file).read().strip()
    if int(seqNum) > 1:
        print('More than one protein sequences found in %s. Please check your protein sequence file.' % protein_file)
        sys.exit(1)

def check_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def check_null_files(infile):
    if os.path.getsize(infile) == 0:
        return True
    else:
        return False

def check_minced_files(indir):
    acceptable = 0
    for minced_result in os.listdir(indir):
        if minced_result.endswith('.minced'):
            filesize = os.path.getsize(os.path.join(indir,minced_result))
            if filesize == 0:
                os.system("rm "+os.path.join(indir,minced_result)+" "+\
                          os.path.join(indir,minced_result.replace('.minced','.minced.gff'))+" "+\
                          os.path.join(indir,minced_result.replace('.minced','.minced.fa')))
            else:
                acceptable += 1
    if acceptable == 0:
        return False
    else:
        return True

def main():

    check_environment()
    check_internet_connection()


if __name__ == '__main__':
    main()
