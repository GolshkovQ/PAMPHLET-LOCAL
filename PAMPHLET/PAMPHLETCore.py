#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
PAMPHLET - PAM Prediction HomoLogous Enhancement Toolkit
'''

import sys
running_python3 = False
if sys.version_info > (3, 0):
    running_python3 = True

import os
import time


print('''

██████╗  █████╗ ███╗   ███╗██████╗ ██╗  ██╗██╗     ███████╗████████╗   ██╗      ██████╗  ██████╗ █████╗ ██╗     
██╔══██╗██╔══██╗████╗ ████║██╔══██╗██║  ██║██║     ██╔════╝╚══██╔══╝   ██║     ██╔═══██╗██╔════╝██╔══██╗██║     
██████╔╝███████║██╔████╔██║██████╔╝███████║██║     █████╗     ██║█████╗██║     ██║   ██║██║     ███████║██║     
██╔═══╝ ██╔══██║██║╚██╔╝██║██╔═══╝ ██╔══██║██║     ██╔══╝     ██║╚════╝██║     ██║   ██║██║     ██╔══██║██║     
██║     ██║  ██║██║ ╚═╝ ██║██║     ██║  ██║███████╗███████╗   ██║      ███████╗╚██████╔╝╚██████╗██║  ██║███████╗
╚═╝     ╚═╝  ╚═╝╚═╝     ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝╚══════╝   ╚═╝      ╚══════╝ ╚═════╝  ╚═════╝╚═╝  ╚═╝╚══════╝                                                                                             

''')

print("PAMPHLET-LOCAL")
print("PAM Prediction HomoLogous Enhancement Toolkit - LOCAL VERSION")

print("===============================================================")
print("Start loading resources...")
resourcesTime = time.time()
runningTime = time.time()

from PAMPHLET import PAMPHLETParams
from PAMPHLET import PAMPHLETReviser
from PAMPHLET import PAMPHLETResources

print("Loading resources finished. Time used: %s seconds" % (time.time() - resourcesTime))

def main():
    print("===============================================================")
    print("Initializing parameters...")

    arg_parser = PAMPHLETParams.get_arg_parser()
    opts = arg_parser.parse_args()

    SpacerFile = opts.spacer
    RepeatSeq = opts.repeat
    ProteinFile = opts.protein
    OutDir = opts.outdir
    FlankLength = opts.flanklen
    CasDatabase = opts.casDB
    BacteriaDatabase = opts.bacteriaDB
    ProkDatabase = opts.protoDB
    BlastMode = opts.blastmode
    ProteinCOV = opts.pcovs
    ProteinIDENT = opts.pident
    RepeatIDENT = opts.rident
    ReviceLenStatus = opts.reviseLen
    FreqMode = opts.freqmode
    orientations = opts.orientation
    BlastThreads = opts.BlastThreads

    ### FILE CHECKS

    PAMPHLETParams.check_spacer_sequences(SpacerFile)
    PAMPHLETParams.check_environment()
    PAMPHLETParams.check_outdir(OutDir)

    ### MODE CHECKS AND SPACER REVISE

    print("Initializing parameters finished. Time used: %s seconds" % (time.time() - runningTime))
    print("===============================================================")

    print("Using protein sequences to revise spacers...")
    PAMPHLETParams.check_protein_files(ProteinFile)
    finalspacer = PAMPHLETReviser.main(SpacerFile,RepeatSeq,ProteinFile,OutDir,BacteriaDatabase,ProteinCOV,ProteinIDENT,RepeatIDENT,CasDatabase)
    print("Using protein sequences to revise spacers finished. Time used: %s seconds" % (time.time() - runningTime))
    print("===============================================================")
    if finalspacer == False:
        finalspacer = SpacerFile

    print("Start running spacer blast...")
    print("This step may take a long time, please wait patiently...")
    print("===============================================================")

    ### RENAME SPACER ID
    SpacerIDLenList = PAMPHLETResources.rename_spacer_id(finalspacer,finalspacer+".label", orientations)
    finalspacer = finalspacer + ".tmp"

    MergedSpacerBlastOutput = os.path.join(OutDir,"MergedSpacerBlastOutput.txt")

    PAMPHLETResources.run_spacer_blast(finalspacer,MergedSpacerBlastOutput,BlastMode,ProkDatabase,BlastThreads)

    ### GET SIGNIFICANT SPACER BLAST OUTPUT AND NON-SIGNIFICANT HIT SPACER ID
    SignificantSpacerBlastOutput = os.path.join(OutDir,"SignificantSpacerBlastOutputRaw.txt")
    InSignificantSpacerID = PAMPHLETResources.get_significant_spacer_blast_output(MergedSpacerBlastOutput,SignificantSpacerBlastOutput, SpacerIDLenList)
    
    ### SET FINAL SPACER HIT FILE
    FinalSpacerHitFile = os.path.join(OutDir,"FinalSpacerHit.txt")

    print("Spacer blast finished. Time used: %s seconds" % (time.time() - runningTime))
    print("===============================================================")

    print("This is the length of InSignificantSpacerID: ", len(InSignificantSpacerID))

    if len(InSignificantSpacerID) != 0 and ReviceLenStatus == True:
        
        print("First round spacer blast finished. Time used: %s seconds" % (time.time() - runningTime))
        print("Start revising spacer length for non-significant spacer blast hits...")
        print("This step may take a long time, please wait patiently...")
        print("Revised spacers: "+"; ".join(InSignificantSpacerID))
        print("===============================================================")

        RevisedLengthSpacerFile = os.path.join(OutDir,"RevisedLengthSpacerFile.txt")
        PAMPHLETResources.revise_spacer_length(finalspacer,RevisedLengthSpacerFile,InSignificantSpacerID)
        RevisedSpacerBlastOutput = os.path.join(OutDir,"RevicedMergedSpacerBlastOutput.txt")
        PAMPHLETResources.run_spacer_blast(RevisedLengthSpacerFile,RevisedSpacerBlastOutput,BlastMode,ProkDatabase)

        RevisedSpacerLen = PAMPHLETResources.build_revised_spacer_len_dict(SpacerIDLenList,InSignificantSpacerID)
        SignificantRevisedSpacerBlastOutput = os.path.join(OutDir,"SignificantRevisedSpacerBlastOutputRaw.txt")
        FalseSpacerIDList = PAMPHLETResources.get_significant_spacer_blast_output(RevisedSpacerBlastOutput,SignificantRevisedSpacerBlastOutput,RevisedSpacerLen)

        print("Revised false spacer: "+"; ".join(FalseSpacerIDList))
        print("Getting final spacer blast output...")

        if os.path.getsize(SignificantRevisedSpacerBlastOutput) == 0:
            print("WARNING: No significant spacer blast hit for revised spacer length.")
            if os.path.getsize(SignificantSpacerBlastOutput) == 0:
                print("No significant spacer blast hit for both original spacer and revised spacer .")
                print("Cannot find spacer PAM sequence.")
                sys.exit(1)
            else:
                print("Use original spacer blast hit.")
                ### Also here I change this line from mv to cp, fuck this code
                os.system("mv %s %s" % (SignificantSpacerBlastOutput,FinalSpacerHitFile))
        else:
            PAMPHLETResources.merge_spacer_blast_output(SignificantSpacerBlastOutput,SignificantRevisedSpacerBlastOutput,FinalSpacerHitFile)

    else:
        print("Use original spacer blast hit.")
        print("Spacer blast finished. Time used: %s seconds" % (time.time() - runningTime))
        print("===============================================================")
        ### haha i change this line from mv to cp
        os.system("cp %s %s" % (SignificantSpacerBlastOutput,FinalSpacerHitFile))

    if os.path.getsize(FinalSpacerHitFile) == 0:
        print("No significant spacer blast hit found. PAMPHLET finished.")
        sys.exit(1)

    ### GET GENOME SEQDUMP FILE
    print("Extract genome sequence based on significant spacer blast output...")
    GenomeSeqdumpFile = os.path.join(OutDir,"GenomeSeqdump.txt")
    GenomeCollection = ProkDatabase+".fasta"
    PAMPHLETResources.get_genome_seqdump_files(FinalSpacerHitFile,GenomeCollection,GenomeSeqdumpFile)
    
    print("Genome seqdump file extraction finished. Time used: %s seconds" % (time.time() - runningTime))
    print("===============================================================")

    ### GET FLANK SEQUENCE BASED ON SIGNIFICANT SPACER BLAST OUTPUT

    print("Extract flank sequence and calculate PAM...")

    UpStreamFlankSeq = os.path.join(OutDir,"UpStreamFlankSeq.txt")
    DownStreamFlankSeq = os.path.join(OutDir,"DownStreamFlankSeq.txt")
    tempMinCEDDir = os.path.join(OutDir,"SpacerSeqdumpMinCED/")
    PAMPHLETParams.check_directory(tempMinCEDDir)
    UpstreamDict, DownstreamDict, SortedUpDupDict, SortedDownDupDict = PAMPHLETResources.get_flank_seq(FinalSpacerHitFile,GenomeSeqdumpFile,UpStreamFlankSeq,DownStreamFlankSeq,FlankLength,tempMinCEDDir)

    if len(UpstreamDict) == 0 or len(DownstreamDict) == 0:
        print("Cannot find upstream or downstream flank sequence.")
        sys.exit(1)

    ### CONVERT SEQUENCE DICT TO FREQUENCY LIST
    UpstreamFreqList = PAMPHLETResources.convert_seq_to_freq(UpstreamDict,FlankLength,FreqMode)
    DownstreamFreqList = PAMPHLETResources.convert_seq_to_freq(DownstreamDict,FlankLength,FreqMode)
    print(UpstreamFreqList)
    print(DownstreamFreqList)

    ### WRITE FREQDICT TO FILE, BUILD MATRIX
    UpstreamFreqFile = os.path.join(OutDir,"UpstreamFreq.txt")
    DownstreamFreqFile = os.path.join(OutDir,"DownstreamFreq.txt")
    ### This is my temp add function to run both positive and negative strand in one run
    NegativeUpstreamFreqFile = os.path.join(OutDir,"NegativeUpstreamFreq.txt")
    NegativeDownstreamFreqFile = os.path.join(OutDir,"NegativeDownstreamFreq.txt")
    PAMPHLETResources.write_freq_dict_to_file(UpstreamFreqList,UpstreamFreqFile)
    PAMPHLETResources.write_freq_dict_to_file(DownstreamFreqList,DownstreamFreqFile)
    PAMPHLETResources.write_freq_dict_to_file_neg(UpstreamFreqList,NegativeDownstreamFreqFile)
    PAMPHLETResources.write_freq_dict_to_file_neg(DownstreamFreqList,NegativeUpstreamFreqFile)

    ### NOW, DRAW WEBLOGO AND SAVE TO FILE
    UpstreamLogoFile = os.path.join(OutDir,"UpstreamLogo.png")
    DownstreamLogoFile = os.path.join(OutDir,"DownstreamLogo.png")
    NegativeUpstreamLogoFile = os.path.join(OutDir,"NegativeUpstreamLogo.png")
    NegativeDownstreamLogoFile = os.path.join(OutDir,"NegativeDownstreamLogo.png")
    PAMPHLETResources.draw_weblogo(UpstreamFreqFile,UpstreamLogoFile,"Upstream",FlankLength)
    PAMPHLETResources.draw_weblogo(DownstreamFreqFile,DownstreamLogoFile,"Downstream",FlankLength)
    PAMPHLETResources.draw_weblogo_neg(NegativeUpstreamFreqFile,NegativeUpstreamLogoFile,"Upstream",FlankLength)
    PAMPHLETResources.draw_weblogo_neg(NegativeDownstreamFreqFile,NegativeDownstreamLogoFile,"Downstream",FlankLength)

    ### PAMPHLET FINISHED
    print("PAMPHLET finished.")
    print("Total time used: %s seconds" % (time.time() - runningTime))
    print("===============================================================")

    if len(SortedUpDupDict) <= 10:
        TopSize = len(SortedUpDupDict)
    else:
        TopSize = 10
    for i in range(0,TopSize):
        print("RANK "+str(i+1)+":\n")
        print("Upstream: "+str(SortedUpDupDict[i][1])+"\t"+str(SortedUpDupDict[i][0]))
        print("Downstream: "+str(SortedDownDupDict[i][1])+"\t"+str(SortedDownDupDict[i][0]))
        print("===============================================================")

if __name__ == "__main__":
    main()
