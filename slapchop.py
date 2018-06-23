#!/usr/bin/env python3

# Simply Looking At Pairwise Comparisons to Optimally Parition
#
# Modifying to be more general.
#

import re
import multiprocessing 
import itertools
import argparse
import subprocess
import sys
import time
import os.path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

def alignChop(inputQueue,lock):
    while True:
        inputEntryList = inputQueue.get()
        if inputEntryList == None:
            print(multiprocessing.current_process().name+\
                str(multiprocessing.current_process().pid)+" is done",file=open(args.logFile,"a"))
            return 1
        else:
            outputPass = []
            outputFail = []
            outputReport = []
            print("\t"+multiprocessing.current_process().name+\
                " chewing bite of "+str(len(inputEntryList))+" entries",file=open(args.logFile,"a"))
            for inputEntry in inputEntryList:
                idline = inputEntry[0].rstrip().split(" ")
                inputRecord=SeqRecord(Seq(inputEntry[1].rstrip()),id=idline[0])
                inputRecord.letter_annotations['phred_quality'] = \
                    inputEntry[3][0:len(inputEntry[1].rstrip())]
        
                choppa = []
            
                aln1 = pairwise2.align.localmd(
                    inputRecord.seq, args.fixed1pattern,
                    args.matchScore,args.mismatchScore,
                    args.readGapOpenScore,args.readGapExtendScore,
                    args.templateGapOpenScore,args.templateGapExtendScore,
                    penalize_end_gaps=(True,False),
                    one_alignment_only=True
                )
            #    print(pairwise2.format_alignment(*aln[0]))
            #    print(aln[0])
                Score1 = round(aln1[0][2],3)
                AlignmentStart1 = aln1[0][3]
                AlignmentEnd1 = aln1[0][4]
            
                aln2 = pairwise2.align.localmd(
                    inputRecord.seq, args.fixed2pattern,
                    args.matchScore,args.mismatchScore,
                    args.readGapOpenScore,args.readGapExtendScore,
                    args.templateGapOpenScore,args.templateGapExtendScore,
                    penalize_end_gaps=(True,False),
                    one_alignment_only=True
                )
            #    print(pairwise2.format_alignment(*aln[0]))
            #    print(aln[0])
                Score2 = round(aln2[0][2],3)
                AlignmentStart2 = aln2[0][3]
                AlignmentEnd2 = aln2[0][4]
            
                indexSeq = inputRecord[0:AlignmentStart1]
                fixed1 = inputRecord[AlignmentStart1:AlignmentEnd1]
                strainSeq = inputRecord[(AlignmentEnd1-3):(AlignmentEnd1+25)]
                fixed2 = inputRecord[AlignmentStart2:AlignmentEnd2]
                tailSeq = inputRecord[AlignmentStart2:]
            
                alnUMI = pairwise2.align.localmd(
                    tailSeq.seq,args.umipattern,
                    args.matchScore,args.mismatchScore,
                    args.readGapOpenScore,args.readGapExtendScore,
                    args.templateGapOpenScore,args.templateGapExtendScore,
                    penalize_end_gaps=(True,False),
                    one_alignment_only=True
                )
            #    print(pairwise2.format_alignment(*aln[0]))
            #    print(aln[0])
                ScoreUMI = round(alnUMI[0][2],3)
                AlignmentStartUMI = alnUMI[0][3]
                AlignmentEndUMI = alnUMI[0][4]
        
                umiSeq = tailSeq[(alnUMI[0][3]):(alnUMI[0][4])]
            #    umiSeq = umiSeq[16:31:1]
            #    umiSeq = umiSeq[17:31:2]
                umiSeq = umiSeq[18:29:2]
            
                outputString = str(strainSeq.id)+"_"+str(umiSeq.seq)+"\n"+\
                    str(indexSeq.seq)+str(strainSeq.seq)+"\n"+"+"+"\n"+\
                    indexSeq.letter_annotations['phred_quality']+\
                    strainSeq.letter_annotations['phred_quality']+"\n"
            
                try:
                    if eval(args.filters):
                        outputPass.append(outputString)
                    else:
                        outputFail.append(str(inputRecord.id)+"\n"+\
                            str(inputRecord.seq)+"\n"+"+"+"\n"+\
                            inputRecord.letter_annotations['phred_quality']+"\n")
                except:
                    outputFail.append(outputString)
            
                outputReport.append(inputRecord.id+"	"+
                    str(Score1)+"	"+str(Score2)+"	"+str(ScoreUMI)+"	"+
                    str(AlignmentStart1)+"	"+str(AlignmentEnd1)+"	"+
                    str(AlignmentStart2)+"	"+str(AlignmentEnd2)+"	"+
                    str(AlignmentStartUMI)+"	"+str(AlignmentEndUMI)+"	"+
                    str(indexSeq.seq)+"	"+
                    str(fixed1.seq)+"	"+
                    str(strainSeq.seq)+"	"+
                    str(fixed2.seq)+"	"+
                    str(tailSeq.seq)+"\n")
            
            lock.acquire()
            print("\t\t"+multiprocessing.current_process().name+\
                " writing out bite of "+str(len(inputEntryList))+" entries",file=open(args.logFile,"a"))
            outputReportHandle = open(args.outputReport,"a")
            for i in outputReport:
                outputReportHandle.write(i)
            outputReportHandle.close()
            if args.filters:
                outputPassHandle = open(args.outputBase+"_pass.fastq","a")
                for i in outputPass:
                    outputPassHandle.write(i)
                outputPassHandle.close()
                outputFailHandle = open(args.outputBase+"_fail.fastq","a")
                for i in outputFail:
                    outputFailHandle.write(i)
                outputFailHandle.close()
            else:
                outputFailHandle = open(args.outputBase+"_all.fastq","a")
                for i in outputFail:
                    outputFailHandle.write(i)
                outputFailHandle.close()
            lock.release()


def reader(queue,inputFastqFileName,biteSize):
    with open(inputFastqFileName,'r') as f:
        print(multiprocessing.current_process().name+" trying to read "+\
            inputFastqFileName+" into memory",file=open(args.logFile,"a"))
#        inputFastqList = f.readlines()
        iterf = iter(f)
        while True:
            i = 0
            thisBite = []
            while i < biteSize:
#                anEntry = inputFastqList[0:4]
#                del inputFastqList[0:4]
#                if anEntry:
                try:
                    anEntry = [next(iterf),next(iterf),next(iterf),next(iterf)]
                    thisBite.append(anEntry)
                    i += 1
                except:
                    break
            if thisBite:
                print("Queue is this big, approximately: "+str(sys.getsizeof(thisBite)*queue.qsize()),file=open(args.logFile,"a"))
                while (sys.getsizeof(thisBite)*queue.qsize()/1000000) > int(args.maxQueueSize):
                    print("Waiting for 10 seconds for queue to get chewed down a bit...",file=open(args.logFile,"a"))
                    time.sleep(10)
                queue.put(thisBite)
                print(multiprocessing.current_process().name+\
                    " took a bite of size "+str(len(thisBite))+", put it on queue( size "+\
                    str(queue.qsize())+" bites)",file=open(args.logFile,"a"))
            else:
                for i in range(int(args.processes)+1):
                    queue.put(None)
                break
    return 1



#####
# Main script
#####

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=""+
        "slapchop.py")
    parser.add_argument("input-fastq",
        help="The FASTQ formatted file to process.")
    parser.add_argument("--processes",default=1)
    parser.add_argument("--bite-size",
        help="The size of bites to chomp off of the input file for multi-process",
        default=1000)
#
    parser.add_argument("--operation","-o",action="append",
        help="The pattern to match and extract. This has a "+
            "specific and sort of complicated syntax. To remind "+
            "you, its pretty much like: ''."+
            "YOU CAN STORE MULTIPLE, and just wire the outputs to "+
            "inputs to build chains of operations.")
    parser.add_argument("--filter","-f",action="append",
        help="A filter for eliminating reads that don't pass some "+
            "alignment based cutoff. This is specified per "+
            "operation, so remember the name from above. Syntax: ''")
#
    parser.add_argument("--match-score",default=1)
    parser.add_argument("--mismatch-score",default=0.001,
        help="Has to be above 0 in order to use local alignments.")
    parser.add_argument("--gap-open",default=-1)
    parser.add_argument("--gap-extend",default=-1)
    parser.add_argument("--read-gap-open",default=None)
    parser.add_argument("--read-gap-extend",default=None)
    parser.add_argument("--regex-gap-open",default=None)
    parser.add_argument("--regex-gap-extend",default=None)
#
    parser.add_argument("output-base",
        help="Base name for the output, will be used to make, "+
            "for example: basename.fastq, basename.report")
    parser.add_argument("--write-report",action='store_true',
        help="Add this flag to print a report of per-read "+
            "statistics, that's a lot of disk writes btw.")
    parser.add_argument("--maxQueueSize",help="in gigs",default=10)
#
    args=parser.parse_args()

    if args.read_gap_open is None: args.read_gap_open = args.gap_open
    if args.read_gap_extend is None: args.read_gap_extend = args.gap_extend
    if args.regex_gap_open is None: args.regex_gap_open = args.gap_open
    if args.regex_gap_extend is None: args.regex_gap_extend = args.gap_extend

#'opname: [INPUT] > ATCG(?P<regex>regularExpression)With(?P<capgroups>NamedCaptureGroups)ATCG > [regex,capgroups]' 
#'opname: '

#[('input-fastq', 'zerp'), ('processes', 1), ('bite_size', 1000), ('operation', ['XYZ', '123']), ('filter', None), ('mismatch_score', 0.001), ('gap_open', -1), ('gap_extend', -1), ('read_gap_open', -1), ('read_gap_extend', -1), ('regex_gap_open', -1), ('regex_gap_extend', -1), ('output-base', 'output'), ('write_report', False), ('maxQueueSize', 10)])
    
    # FIRST, report what we're doing.
    print()
    print("I'm reading in '"+vars(args)["input-fastq"]+"', "+
        "applying these operations of alignment:")
    try:
        for i in args.operation:
            print("\t"+i)
    except:
        print("Wait a second, there's no operations to be done! "+
            "Exiting...")
        exit(1)
    print("... with these filters:")
    try:
        for i in args.filter:
            print("\t"+i)
    except:
        print("\t( # no filters defined )")
    print()
    print("I'm going to use these parameters for the alignments:")
    print("Match "+str(args.match_score)+", "+
        "mismatch "+str(args.mismatch_score)+", "+
        "read gap open "+str(args.read_gap_open)+", "+
        "read gap extend "+str(args.read_gap_extend)+", "+
        "regex gap open "+str(args.regex_gap_open)+", "+
        "regex gap extend "+str(args.regex_gap_extend)+".")
    print()
    print("Then, I'm going to write out a FASTQ file to '"+
        vars(args)["output-base"]+".fastq'",end="")
    if args.write_report:
        print(" and a report to '"+
            vars(args)["output-base"]+"_report.csv'",end="")
    print(".")

    print()
    print("I will proceed to process the file with "+
        str(args.processes)+" processes operating in chunks of "+
        str(args.bite_size)+" records.")
    
    print()
    print("BEGIN")
    print()
    
    if os.path.isfile(vars(args)["output-base"]+".fastq"):
        print("File "+vars(args)["output-base"]+".fastq "+
            "exits, so I'm quitting before you ask me to do "+
            "something you'll regret.")
    if os.path.isfile(vars(args)["output-base"]+"_report.csv"):
        print("File "+vars(args)["output-base"]+".fastq "+
            "exits, so I'm quitting before you ask me to do "+
            "something you'll regret.")

    biteSize = int(args.biteSize)

    manager = multiprocessing.Manager()
    queue = manager.Queue()
    lock = manager.Lock()

    jobs = []
    jobs.append(multiprocessing.Process(target=reader,\
            args=(queue,args.inputFastq,biteSize),name="TheReader"))
    for i in range(1,int(args.processes)):
        jobs.append(multiprocessing.Process(target=alignChop,args=(queue,lock),\
            name="TheChewer"+str(i)))

    for i in jobs:
        i.start()
        print(i,file=open(args.logFile,"a"))

    jobs[0].join()
    print("TheReader is done, allocating another Chewer number 0",file=open(args.logFile,"a"))
    jobs[0] = multiprocessing.Process(target=alignChop,args=(queue,lock),\
        name="TheChewer0")
    jobs[0].start()

    for i in jobs:
        if i.is_alive():
            i.join()

    print("We done here?",file=open(args.logFile,"a"))
    for i in jobs:
        print(i,file=open(args.logFile,"a"))
