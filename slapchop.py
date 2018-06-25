#!/usr/bin/env python3

# Simply Looking At Pairwise Comparisons to Optimally Parition
#
# Modifying to be more general.
#
#usage for testing
# rm sobaseqout*; ./slapchop.py example_sobaseq.fastq sobaseqout --operation "get_sample: INPUT > (?P<sample>[ATCG]{5})(GTCCTCGAGGTCTCT){e<=2}(?P<rest>.*)$" -o "get_strain: rest > (?P<strain>[ATCG]{10,26})(GCGTACGCTGCAGGT){e<=2}.*" --filter "sample_length == 5" --write-report --bite-size 10 --processes 2; wc -l sobaseqout*

import re
import multiprocessing 
import itertools
import argparse
import subprocess
import sys
import time
import os.path
from Bio import Seq, SeqRecord, pairwise2, SeqIO
import regex

def reader(input_line_queue,input_fastq,
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
    seq_holder['INPUT'] = inputRecord 
#rewrite as a class

    for operation_name, operation in operations_dict.items():

        try:
            if len(seq_holder[operation[0]]) == 0:
                continue
        except:
            continue

        fuzzy_match = regex.search(
            operation[1], # the seq_pattern to match
            str(seq_holder[operation[0]].seq), # the input seq 
            regex.BESTMATCH )

        if fuzzy_match is None:
            continue
        else:
            (scores_holder[operation_name+'_substitutions'],
                scores_holder[operation_name+'_insertions'],
                scores_holder[operation_name+'_deletions']
                ) = fuzzy_match.fuzzy_counts
            for match_name in fuzzy_match.groupdict():
                seq_holder[match_name] = \
                    seq_holder[operation[0]]\
                    [slice(*fuzzy_match.span(match_name))]
                (scores_holder[match_name+'_start'],
                    scores_holder[match_name+'_end']
                    ) = fuzzy_match.span(match_name)
                scores_holder[match_name+'_length'] = \
                    (scores_holder[match_name+'_end'] - 
                        scores_holder[match_name+'_start'])

    evaluated_filters = evaluate_filters(args.filter,scores_holder)

    if not all(evaluated_filters):
        pass
    else:
        print("WRITE OUT")
#        outputString = str(strainSeq.id)+"_"+str(umiSeq.seq)+"\n"+\
#            str(indexSeq.seq)+str(strainSeq.seq)+"\n"+"+"+"\n"+\
#            indexSeq.letter_annotations['phred_quality']+\
#            strainSeq.letter_annotations['phred_quality']+"\n"
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

def evaluate_filters(filters,scores_holder):
    locals().update(scores_holder)
    return_object = []
    try:
        for each_filter in filters:
            if eval(each_filter):
                return_object.append(True)
            else:
                return_object.append(each_filter)
    except:
        return([False])
    return(return_object)





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
    parser.add_argument("--output-format",
        help="format for the output file")
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

#####

    operations_dict = dict()

    for each in args.operation:

        (name, instruction) = each.split(":")
        (input_string, regex_string) = instruction.strip().split(" > ")
        operations_dict[name] = [input_string, regex_string]

    # FIRST, report what we're doing.
    print()
    print("I'm reading in '"+vars(args)["input-fastq"]+"', "+
        "applying these operations of alignment:")
    try:
        for key, value in operations_dict.items():
            print("\t"+key+":"+
                "\n\t\tfrom : "+value[0]+
                "\n\t\textract groups with regex : '"+value[1]+"'"+
                "\n")
    except:
        print("Wait a second, there's no operations to be done! "+
            "Exiting...")
        exit(1)

#####

    print("...and with these filters:")
    try:
        for i in args.filter:
            print("\t"+i)
    except:
        print("\t( # no filters defined )")

#####

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
        jobs.append(multiprocessing.Process(target=reader,
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
