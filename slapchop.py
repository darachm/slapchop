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
from Bio import Seq, SeqRecord, pairwise2, SeqIO

def alignChopRead(input_line_queue,input_fastq,
                    out_lock,output_fastq,
                    report_lock,report_csv,
                    bite_size,operations_dict):

    while True:

        current_pos = input_line_queue.get()

        if current_pos == "poisonpill":
            input_line_queue.put("poisonpill")
            return(0)

        report_lock.acquire()
        print(multiprocessing.current_process().name+
            " trying to read into memory chunk of size "+
            str(bite_size)+" records.",
            file=open(report_csv,"a"))
        report_lock.release()

        ifqp = open(input_fastq,"r")
        ifqp.seek(current_pos)
        chunk = []
        for i in range(bite_size*4):
            chunk.append(ifqp.readline())
#                try:
#                    anEntry = [next(iterf),next(iterf),next(iterf),next(iterf)]
#                    thisBite.append(anEntry)
#                    i += 1
#                except:
#                    break


        current_pos = ifqp.tell()
        if ifqp.readline() == "":
            input_line_queue.put("poisonpill")
            report_lock.acquire()
            print(multiprocessing.current_process().name+\
                str(multiprocessing.current_process().pid)+
                " is done",file=open(report_csv,"a"))
            report_lock.release()
        else:
            input_line_queue.put(current_pos)

        for slice_base in itertools.islice(range(len(chunk)),0,None,4):

            if chunk[slice_base] == "":
                break
    
            alignChop(chunk[slice(slice_base,(slice_base+4))],
                operations_dict)

            out_lock.acquire()
#                print(i.strip(),file=open(output_fastq,"a"))
            out_lock.release()


def alignChop(record,operations_dict):
    inputRecord = SeqRecord.SeqRecord(Seq.Seq(record[1].rstrip()),
        id = record[0].rstrip().split(" ")[0])
    inputRecord.letter_annotations['phred_quality'] = \
        record[3][0:len(record[1].rstrip())]

    scores_holder = dict()
    seq_holder = dict()
    seq_holder['INPUT'] = [ inputRecord ]
#rewrite as a class

    print()

    for operation_name, operation in operations_dict.items():

        try:
            if len(seq_holder[operation[0].upper()]) == 0:
                continue
        except:
            continue

        aln = pairwise2.align.localmd(
            seq_holder[operation[0].upper()][0].seq, # the input seq
            operation[1], # the seq_pattern to match
            args.match_score, args.mismatch_score,
            args.read_gap_open, args.read_gap_extend,
            args.seq_pattern_gap_open, args.seq_pattern_gap_extend,
            penalize_end_gaps=(True,False),
            one_alignment_only=True
        )[0]

        scores_holder[operation_name] = aln[2:len(aln)]

        regex_match = operation[2].search(aln[0][aln[3]:aln[4]])
        if regex_match is None:
            continue
        else:
            for match_name in regex_match.groupdict():
                seq_holder[match_name] = \
                    seq_holder[operation[0].upper()][0]\
                    [slice(*regex_match.span(match_name))]
    

    print(seq_holder)


#            #    print(pairwise2.format_alignment(*aln[0]))
#            #    print(aln[0])
#                Score1 = round(aln1[0][2],3)
#                AlignmentStart1 = aln1[0][3]
#                AlignmentEnd1 = aln1[0][4]
#            
#                aln2 = pairwise2.align.localmd(
#                    inputRecord.seq, args.fixed2pattern,
#                    args.matchScore,args.mismatchScore,
#                    args.readGapOpenScore,args.readGapExtendScore,
#                    args.templateGapOpenScore,args.templateGapExtendScore,
#                    penalize_end_gaps=(True,False),
#                    one_alignment_only=True
#                )
#            #    print(pairwise2.format_alignment(*aln[0]))
#            #    print(aln[0])
#                Score2 = round(aln2[0][2],3)
#                AlignmentStart2 = aln2[0][3]
#                AlignmentEnd2 = aln2[0][4]
#            
#                indexSeq = inputRecord[0:AlignmentStart1]
#                fixed1 = inputRecord[AlignmentStart1:AlignmentEnd1]
#                strainSeq = inputRecord[(AlignmentEnd1-3):(AlignmentEnd1+25)]
#                fixed2 = inputRecord[AlignmentStart2:AlignmentEnd2]
#                tailSeq = inputRecord[AlignmentStart2:]
#            
#                alnUMI = pairwise2.align.localmd(
#                    tailSeq.seq,args.umipattern,
#                    args.matchScore,args.mismatchScore,
#                    args.readGapOpenScore,args.readGapExtendScore,
#                    args.templateGapOpenScore,args.templateGapExtendScore,
#                    penalize_end_gaps=(True,False),
#                    one_alignment_only=True
#                )
#            #    print(pairwise2.format_alignment(*aln[0]))
#            #    print(aln[0])
#                ScoreUMI = round(alnUMI[0][2],3)
#                AlignmentStartUMI = alnUMI[0][3]
#                AlignmentEndUMI = alnUMI[0][4]
#        
#                umiSeq = tailSeq[(alnUMI[0][3]):(alnUMI[0][4])]
#            #    umiSeq = umiSeq[16:31:1]
#            #    umiSeq = umiSeq[17:31:2]
#                umiSeq = umiSeq[18:29:2]
#            
#                outputString = str(strainSeq.id)+"_"+str(umiSeq.seq)+"\n"+\
#                    str(indexSeq.seq)+str(strainSeq.seq)+"\n"+"+"+"\n"+\
#                    indexSeq.letter_annotations['phred_quality']+\
#                    strainSeq.letter_annotations['phred_quality']+"\n"
#            
#                try:
#                    if eval(args.filters):
#                        outputPass.append(outputString)
#                    else:
#                        outputFail.append(str(inputRecord.id)+"\n"+\
#                            str(inputRecord.seq)+"\n"+"+"+"\n"+\
#                            inputRecord.letter_annotations['phred_quality']+"\n")
#                except:
#                    outputFail.append(outputString)
#            
#                outputReport.append(inputRecord.id+"	"+
#                    str(Score1)+"	"+str(Score2)+"	"+str(ScoreUMI)+"	"+
#                    str(AlignmentStart1)+"	"+str(AlignmentEnd1)+"	"+
#                    str(AlignmentStart2)+"	"+str(AlignmentEnd2)+"	"+
#                    str(AlignmentStartUMI)+"	"+str(AlignmentEndUMI)+"	"+
#                    str(indexSeq.seq)+"	"+
#                    str(fixed1.seq)+"	"+
#                    str(strainSeq.seq)+"	"+
#                    str(fixed2.seq)+"	"+
#                    str(tailSeq.seq)+"\n")
#            
#            lock.acquire()
#            print("\t\t"+multiprocessing.current_process().name+\
#                " writing out bite of "+str(len(inputEntryList))+" entries",file=open(args.logFile,"a"))
#            outputReportHandle = open(args.outputReport,"a")
#            for i in outputReport:
#                outputReportHandle.write(i)
#            outputReportHandle.close()
#            if args.filters:
#                outputPassHandle = open(args.outputBase+"_pass.fastq","a")
#                for i in outputPass:
#                    outputPassHandle.write(i)
#                outputPassHandle.close()
#                outputFailHandle = open(args.outputBase+"_fail.fastq","a")
#                for i in outputFail:
#                    outputFailHandle.write(i)
#                outputFailHandle.close()
#            else:
#                outputFailHandle = open(args.outputBase+"_all.fastq","a")
#                for i in outputFail:
#                    outputFailHandle.write(i)
#                outputFailHandle.close()
#            lock.release()



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
    parser.add_argument("--seq_pattern-gap-open",default=None)
    parser.add_argument("--seq_pattern-gap-extend",default=None)
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
    if args.seq_pattern_gap_open is None: args.seq_pattern_gap_open = args.gap_open
    if args.seq_pattern_gap_extend is None: args.seq_pattern_gap_extend = args.gap_extend

    args.bite_size = int(args.bite_size)

    operations_dict = dict()

    for each in args.operation:

        (name, instruction) = each.split(":")
        (input_string, seq_pattern) = instruction.strip().split(" > ")
        converted_to_seq = re.sub(
            "[^ATCGN]",
            "",
            re.sub(
                "\(\?P<[^>]+>([^\)]+)\)",
                "\\1",
                seq_pattern.upper())
            ) 
        converted_to_regex = seq_pattern.upper()
        while re.search("N[^>]*(?:[^>]|$)",converted_to_regex) is not None:
            converted_to_regex = re.sub("N([^>]*(?:[^>]|$))","[ACTG]\\1",converted_to_regex)
        operations_dict[name] = [input_string, 
            converted_to_seq,
            re.compile(converted_to_regex)]

    # FIRST, report what we're doing.
    print()
    print("I'm reading in '"+vars(args)["input-fastq"]+"', "+
        "applying these operations of alignment:")
    try:
        for key, value in operations_dict.items():
            print("\t"+key+":"+
                "\n\t\tfrom : "+value[0]+
                "\n\t\tsequence pattern : "+value[1]+
                "\n\t\textract groups with regex : "+str(value[2])+
                "\n"
                )
    except:
        print("Wait a second, there's no operations to be done! "+
            "Exiting...")
        exit(1)
    print("...and with these filters:")
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
        "pattern gap open "+str(args.seq_pattern_gap_open)+", "+
        "pattern gap extend "+str(args.seq_pattern_gap_extend)+".")
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
        exit(1)
    if os.path.isfile(vars(args)["output-base"]+"_report.csv"):
        print("File "+vars(args)["output-base"]+".fastq "+
            "exits, so I'm quitting before you ask me to do "+
            "something you'll regret.")
        exit(1)

    #####
    # Multi proc
    #####

    manager  = multiprocessing.Manager()
    input_line_queue = multiprocessing.Queue()
    input_line_queue.put(0)
    (out_lock, report_lock) = (manager.Lock(), manager.Lock())

    jobs = []
    for i in range(1,int(args.processes)+1):
        jobs.append(multiprocessing.Process(target=alignChopRead,
            args=(input_line_queue,vars(args)["input-fastq"],
                out_lock,vars(args)["output-base"]+".fastq",
                report_lock,vars(args)["output-base"]+"_report.csv",
                args.bite_size,operations_dict),
            name="Comrade"+str(i)))
        print("Comrade"+str(i))

    for i in jobs:
        i.start()

    for i in jobs:
        if i.is_alive():
            i.join()

    print("We done here?")
