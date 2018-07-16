#!/usr/bin/env python3

# 
# description to be written ...
#

import re
import multiprocessing 
import itertools
import argparse
#import subprocess
#import sys
#import time
import os.path
from Bio import Seq, SeqRecord, SeqIO
import regex
import json
import collections

def reader(input_line_queue,input_fastq,
                    pass_lock,pass_fastq,
                    fail_lock,fail_fastq,
                    report_lock,report_txt,
                    record_lock,record_csv,
                    bite_size,operations_dict):

    while True:

        current_pos = input_line_queue.get()

        try:
            if current_pos == "poisonpill":
                report_lock.acquire()
                print(multiprocessing.current_process().name+
                        " with pid "+
                        str(multiprocessing.current_process().pid)+
                        " found a poison pill and is done",
                    file=open(report_txt,"a"))
                report_lock.release()
                input_line_queue.put("poisonpill")
                return(0)
        except:
            try:
                report_lock.release()
            except:
                return(0)
            return(0)

        try:
            report_lock.acquire()
            print(multiprocessing.current_process().name+
                " trying to read a chunk of size "+
                str(bite_size)+" fastq records.",
                file=open(report_txt,"a"))
            report_lock.release()
        except:
            try:
                report_lock.release()
            except:
                pass
            pass

        ifqp = open(input_fastq,"r")
        ifqp.seek(current_pos)
        chunk = []
        for i in range(bite_size*4):
            chunk.append(ifqp.readline())

        current_pos = ifqp.tell()
        if ifqp.readline() == "":
            input_line_queue.put("poisonpill")
            report_lock.acquire()
            print(multiprocessing.current_process().name+
                    " with pid "+
                    str(multiprocessing.current_process().pid)+
                    " found the end and is done",
                file=open(report_txt,"a"))
            report_lock.release()
            return(0)
        else:
            input_line_queue.put(current_pos)

        pass_records = []
        fail_records = []
        report_records = []

        for slice_base in itertools.islice(range(len(chunk)),0,None,4):

            if chunk[slice_base] == "":
                break
    
            (passed, output_record, report_object) = \
                alignChop(chunk[slice(slice_base,(slice_base+4))],
                    operations_dict,report_txt)

            if passed:
                pass_records.append(output_record)
            else:
                fail_records.append(output_record)

            report_records.append(report_object)

        
        pass_lock.acquire()
        with open(pass_fastq,"a") as f:
            for i in pass_records:
                print(str(i.id)+"\n"+str(i.seq)+"\n"+"+"+"\n"+
                        i.letter_annotations['phred_quality']+"\n",
                    file=f)
        pass_lock.release()

        fail_lock.acquire()
        with open(fail_fastq,"a") as f:
            for i in fail_records:
                print(str(i.id)+"\n"+str(i.seq)+"\n"+"+"+"\n"+
                        i.letter_annotations['phred_quality']+"\n",
                    file=f)
        fail_lock.release()

        record_lock.acquire()
        with open(record_csv,"a") as f:
            for i in report_records:
                print(i,file=f)
        record_lock.release()

        if args.debug:
            return(0)


def alignChop(record,operations_dict,report_txt):

    input_record = SeqRecord.SeqRecord(Seq.Seq(record[1].rstrip()),
        id = record[0].rstrip().split(" ")[0])
    input_record.letter_annotations['phred_quality'] = \
        record[3][0:len(record[1].rstrip())]

    scores_holder = dict()
    seq_holder = dict()
    seq_holder['input'] = input_record 

#rewrite as a class ????

    try:
        if args.debug:
            report_lock.acquire()
            print(multiprocessing.current_process().name+
                " processing:\n"+
                input_record.seq,
                file=open(report_txt,"a"))
            report_lock.release()
    except:
        try:
            report_lock.release()
        except:
            pass
        pass

    for operation_name, operation in operations_dict.items():

        try:
            if args.debug:
                report_lock.acquire()
                print(multiprocessing.current_process().name+
                    " attempting to match :\n"+
                    operation[1],end="",
                    file=open(report_txt,"a"))
                print(" against "+
                    str(seq_holder[operation[0]].seq),
                    file=open(report_txt,"a"))
                report_lock.release()

            if len(seq_holder[operation[0]]) == 0:
                continue

        except:
            try:
                report_lock.release()
            except:
                pass
            continue

        fuzzy_match = regex.search(
            operation[1], # the seq_pattern to match
            str(seq_holder[operation[0]].seq), # the input seq 
            regex.BESTMATCH )

        try:
            if args.debug:
                report_lock.acquire()
                print(multiprocessing.current_process().name+
                    " match is :\n"+
                    str(fuzzy_match),
                    file=open(report_txt,"a"))
                report_lock.release()
        except:
            try:
                report_lock.release()
            except:
                pass

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

    try:
        if args.debug:
            report_lock.acquire()
            print(multiprocessing.current_process().name+
                " evaluated filters is :\n"+
                evaluated_filters,
                file=open(report_txt,"a"))
            report_lock.release()
    except:
        try:
            report_lock.release()
        except:
            pass
        pass


    if not all(evaluated_filters):
        output_record = input_record
        return((False,output_record,
            "\"FailedFilter\","+
            "\""+input_record.id+"\",\""+
                input_record.seq+"\",\""+
                re.sub("\"","\\\"",re.sub(",","\,",json.dumps(scores_holder)))+"\""
            ))
    else:
        try:
            output_record = evaluate_output_directives(
                args.output_seq,args.output_id,seq_holder) 
            return((True,output_record,
                "\"Passed\","+
                "\""+output_record.id+"\",\""+
                    output_record.seq+"\",\""+
                    re.sub("\"","\\\"",re.sub(",","\,",json.dumps(scores_holder)))+"\""
                ))
        except:
            output_record = input_record
            return((False,output_record,
                "\"FailedEvalOutDirectives\","+
                "\""+input_record.id+"\",\""+
                    input_record.seq+"\",\""+
                    re.sub("\"","\\\"",re.sub(",","\,",json.dumps(scores_holder)))+"\""
                ))


def evaluate_output_directives(output_seq, output_id, seq_holder):
    locals().update(seq_holder)
    return_record = eval(output_seq)
    return_record.id = eval(output_id)
    return(return_record)


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
    parser.add_argument("--debug",action="store_true")
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
    parser.add_argument("--output-id",
        help="format for the output file id, per read that passes "+
            "filter",
        default="input.id")
    parser.add_argument("--output-seq",
        help="format for the output file seq, per read that "+
            "passes filter",
        default="input")
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

#####

    args.bite_size = int(args.bite_size)

#####

    operations_dict = collections.OrderedDict()

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
    (pass_lock, fail_lock, report_lock, record_lock) = (
        manager.Lock(), manager.Lock(), 
        manager.Lock(), manager.Lock() )

    jobs = []
    for i in range(1,int(args.processes)+1):
        jobs.append(multiprocessing.Process(target=reader,
            args=(input_line_queue,vars(args)["input-fastq"],
                pass_lock,vars(args)["output-base"]+"_pass.fastq",
                fail_lock,vars(args)["output-base"]+"_fail.fastq",
                report_lock,vars(args)["output-base"]+"_log.txt",
                record_lock,vars(args)["output-base"]+"_record.csv",
                args.bite_size,operations_dict),
            name="Comrade"+str(i)))
        print("Comrade"+str(i))

    for i in jobs:
        i.start()

    for i in jobs:
        if i.is_alive():
            i.join()

    print("We done here?")
