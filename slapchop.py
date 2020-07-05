#!/usr/bin/env python3

# 
# SLAPCHOP - Simply Looking At Pairwise Comparisons to Help Optimize Parsing
# Name comes from first version that used SW alignments.
# It now uses fuzzy regular expressions, so maybe it should be called something
# else. But who cares.
#

import re
import json
import regex
import os.path
import time
import statistics
import argparse
import itertools
import multiprocessing 
from Bio import Seq, SeqRecord
import tracemalloc
import gzip

import sys
from Bio import SeqIO

def reader(
    input_file,is_gzipped,output_file,failed_file,report_file,
    verbosity,
    #memory_tracking_level,
    operations_array, 
    filters ,
    outputs_array
    ):

    print(
        input_file,is_gzipped,output_file,failed_file,report_file,
        verbosity,
        #memory_tracking_level,
        operations_array, 
        filters ,
        outputs_array
        )

    if input_file is not None:
        if is_gzipped:
            with gzip.open(input_file,"rt") as input_file_gz:
                input_seqs = SeqIO.parse(input_file_gz,"fastq")
        else:
            input_seqs = SeqIO.parse(input_file,"fastq")
    else:
        input_seqs = SeqIO.parse(sys.stdin,"fastq")

    for each_seq in input_seqs:
        chop(each_seq
            )

        (passed, output_record, report_object) = \
            chop(
                chunk[slice(slice_base,(slice_base+4))],
                operations_array, filters,
                output_array
                report, verbosity, memory_tracking_level,
                input_fastq)

        # This is a per-record test
        if passed:
            pass_records.append(output_record)
        else:
            fail_records.append(output_record)
        if if_write_report:
            report_records.append(report_object)

        del passed
        del output_record
        del report_object
    
    # Then we write out the records with the appropriate locks
    with pass_lock:
        with open(pass_fastq,"a") as f:
            for i in pass_records:
                print(str(i.id)+"\n"+str(i.seq)+"\n"+"+"+"\n"+
                        i.letter_annotations['phred_letters'],
                    file=f)

    with fail_lock:
        with open(fail_fastq,"a") as f:
            for i in fail_records:
                print(str(i.id)+"\n"+str(i.seq)+"\n"+"+"+"\n"+
                        i.letter_annotations['phred_letters'],
                    file=f)

    if if_write_report:
        with report_lock:
            with open(report_csv,"a") as f:
                for i in report_records:
                    print(i,file=f)

    # `memory_tracking` memory profiler
    print_memory_tracking(1,"At end of bite",memory_tracking_level)

    return(0)

def chop(
    record, operations_array, filters,
    output_seq_spec, output_id_spec,
    if_write_report, verbosity, memory_tracking_level, input_fastq
    ):

    # `memory_tracking` memory profiler
    print_memory_tracking(2,"At begin of chop()",memory_tracking_level)

    # Making the input record from the raw strings
    r0split = record[0].rstrip().split(" ",maxsplit=1)
    input_record = SeqRecord.SeqRecord(Seq.Seq(record[1].rstrip()),
        id = r0split[0], 
        description = r0split[1] if len(r0split) == 2 else ""
        )
    input_record.letter_annotations['phred_letters'] =  record[3][0:len(record[1].rstrip())]

    # Making a spacer thing
    spacer = SeqRecord.SeqRecord(Seq.Seq("X"),id="spacer")
    spacer.letter_annotations['phred_letters'] = "I"

    # Make the right quality scores
    input_record.letter_annotations['phred_quality'] = [ ord(i)-33 for i in input_record.letter_annotations['phred_letters'] ]

    # We make some holders for these operations
    scores_holder = dict()
    seq_holder = {'spacer': spacer, 'input': input_record}

#rewrite as a class ????

    # Chop grained verbosity
    if verbosity > 2:
        print("\n"+"["+str(time.time())+"]"+" : "+multiprocessing.current_process().name+
            " starting to process :\n  "+
            input_record.id+"\n  "+
            input_record.seq+"\n  "+
            input_record.letter_annotations['phred_letters']
            )

    for each_operation in operations_array:

        # The first element is the name, the next two are used later
        operation_name = each_operation[0]
        operation = each_operation[1:]

        # This should fail if you didn't specify anything taking 
        # from input stream!
        if operation[0] not in seq_holder.keys():
            if verbosity > 3:
                print("\n"+"["+str(time.time())+"]"+" : "+multiprocessing.current_process().name+
                " can't find the sequence named `"+
                operation[0]+"` in the holder, so continuing."
                )
            continue

        if verbosity > 3:
            print("\n"+"["+str(time.time())+"]"+" : "+multiprocessing.current_process().name+
                " attempting to match : "+operation[1]+
                " against "+seq_holder[operation[0]].seq
                )

        # Here we execute the actual meat of the business
        sequence_to_search = str(seq_holder[operation[0]].seq).upper()

        compiled_regex = regex.compile(
            # We use this regex
            operation[1], 
            # And we use the BESTMATCH strategy, I think
            regex.BESTMATCH
            )
        fuzzy_match = compiled_regex.search(
            # to search on this sequence
            sequence_to_search
            )

        if verbosity > 3:
            print("\n"+"["+str(time.time())+"]"+" : "+multiprocessing.current_process().name+
                " match is : "+
                str(fuzzy_match)
                )

        # This is fine, just means the pattern couldn't match at all
        if fuzzy_match is None:
            continue
        # If we did match, then we store them in places
        else:
            # We use tuples to store all the details of the kinds
            # of errors that allowed the match
            (scores_holder[operation_name+'_substitutions'],
                scores_holder[operation_name+'_insertions'],
                scores_holder[operation_name+'_deletions']
                ) = fuzzy_match.fuzzy_counts
            # Then for each of the groups matched by the regex
            for match_name in fuzzy_match.groupdict():
                # We stick into the holder
                # a slice of the input seq, that is the matched
                # span of this matching group
                seq_holder[match_name] = \
                    seq_holder[operation[0]]\
                    [slice(*fuzzy_match.span(match_name))]
                # Then we record the start, end, and length of the
                # matched span
                (scores_holder[match_name+'_start'],
                    scores_holder[match_name+'_end']
                    ) = fuzzy_match.span(match_name)
                scores_holder[match_name+'_length'] = \
                    (scores_holder[match_name+'_end'] - 
                        scores_holder[match_name+'_start'])

    # All these values allow use to apply filters, using this
    # function
    evaluated_filters = evaluate_filters(filters, {**scores_holder, **seq_holder} )

    # This evaluated_filters should be logical list
    if not all(evaluated_filters):

        if verbosity > 2:
            print("\n"+"["+str(time.time())+"]"+" : "+multiprocessing.current_process().name+
                " evaluated the filters as : "+
                str(evaluated_filters)+
                " and so failed."
                )

        # So if we should write this per-record report
        if if_write_report:
            return((False,input_record,
                "\"FailedFilterOnThisInput\","+
                "\""+input_record.id+"\",\""+
                    input_record.seq+"\",\""+
                    re.sub("\"","\\\"",re.sub(",","\,",
                        json.dumps(scores_holder)))+"\""
                ))
        else:
            return((False,input_record,""))
        # If this json dump is empty, it might be because it didn't
        # ever match the first operation, so then just died without
        # building that object

    else:

        try:
            # We attempt to form the correct output record based on
            # the arguments given
            output_record = evaluate_output_directives(
                output_seq_spec,output_id_spec,seq_holder) 

            if verbosity > 2:
                print("\n"+"["+str(time.time())+"]"+" : "+multiprocessing.current_process().name+
                    " evaluated the filters as : "+
                    str(evaluated_filters)+
                    " and so passed!"
                    )

            # If we want to write the report, we make it
            if if_write_report:

                # `memory_tracking` memory profiler
                print_memory_tracking(3,"At end of chop(), while writing",memory_tracking_level)

                return((True,output_record,
                    "\"Passed\","+
                    "\""+output_record.id+"\",\""+
                        output_record.seq+"\",\""+
                        re.sub("\"","\\\"",re.sub(",","\,",
                            json.dumps(scores_holder)))+"\""
                    ))
            # Otherwise, just the record and if it passed
            else:
                return((True,output_record,""))

        except:

            if verbosity > 2:
                print("\n"+"["+str(time.time())+"]"+" : "+multiprocessing.current_process().name+
                    " failed upon forming the output."
                    )

            if if_write_report:
                return((False,input_record,
                    "\"FailedDirectivesToMakeOutputSeq\","+
                    "\""+input_record.id+"\",\""+
                        input_record.seq+"\",\""+
                        re.sub("\"","\\\"",re.sub(",","\,",
                            json.dumps(scores_holder)))+"\""
                    ))
            else:
                return((False,input_record,""))


def evaluate_output_directives(output_seq, output_id, seq_holder):
    # Here we evaluate them but using that dictionary as the global
    # dictionary, because done is better than dogma.
    return_record    = eval(output_seq,{},seq_holder)
    return_record.id = eval(output_id, {},seq_holder)
    return(return_record)


def evaluate_filters(filters,holder):
    return_object = []
    try:
        for each_filter in filters:
            # Here we evaluate them but using that dictionary as the
            # global dictionary, because done is better than dogma.
            if eval(each_filter,globals(),holder):
                return_object.append(True)
            else:
                return([False])
    except:
        return([False])
    return(return_object)


if __name__ == '__main__':

#### defining arguments with argparse module

    # Name and description
    parser = argparse.ArgumentParser(description=""+
        "slapchop.py - we expect to have a standard FASTQ piped in on STDIN"+
        "")

    # Optional input file, but we're expecting STDIN as inputs
    parser.add_argument("--input",default=None,
        help="An FASTQ(Z) file to parse, we expect STDIN but it's an option.")
    # gzip format
    parser.add_argument("-z","--gzipped",action="store_true",
        help="Is input file actually a FASTQZ, so gzipped FASTQ?")

    # Output files
    parser.add_argument("--output",default=None,
        help="Name of output file, incase you want to write it to a file, "+
            "but we expect it to go to STDOUT. "+
            "For example: basename.fastq, basename.report")
    parser.add_argument("--failed",default=None,
        help="Name of output file for failed reads, for debugging."+
            "If you say 'STDOUT' then it'll go there. "
        )
    parser.add_argument("--report",default=None,
        help="Add this flag with a filename to print a report of per-read "+
            "statistics, that's a lot of disk writes btw, but "+
            "would be good in combination with a limited argument "+
            "so that you can spec out the kinds of noise you got "+
            "in your data.")

    # verbosity
    parser.add_argument("-v","--verbose",action="count",default=0,
        help="How much debugging issues should I pipe out to STDERR?"+
            "0 is nothing, "+
            "1 is setup messages and start-stop messsages, "+
            "2 is worker-level details, "+
            "3 is chop-level details, "+
            "4 is each operation level details."
            )

    # Operations
    parser.add_argument("--operation","-o",action="append",
        help="The pattern to match and extract. This has a "+
            "specific and sort of complicated syntax. Refer to "+
            "the documentation via the README.md file, or to "+
            "the original regex module documentation. "+
            "They are chained in the order you specify. "+
            "You can't name a group 'input' or 'spacer' !")

    # Filter specification
    parser.add_argument("--filter","-f",action="append",
        help="Filters for eliminating reads that don't pass some "+
            "alignment based cutoff. This is specified per "+
            "operation, so remember the name from above. "+
            "Every sequence has a _start, _end, and _length "+
            "calculated. You're welcome."
            )

    # Outputs
    parser.add_argument("--output-id",action="append",
        help="A list of output IDs, in the same order as for output-seq. "+
            "This is for reads that pass filter. This is evaluated, so you "+
            "can access .id or .seq properties of the sequences. You "+
            "probably want to include `input.id` in there.",
        )
    parser.add_argument("--output-seq",action="append",
        help="A list of output seqs, in the same order as for output-ids. "+
            "This is a Biopython SeqRecord, so "+
            "you need to just specify the names of the captured "+
            "groups. You can't access .id or _length properties or "+
            "the like! This is not for that, put it in the "+
            "output-id. Why? Because need to access qualities for printing.",
        )

#### Parse, clean up, and possibly report arguments

    args = parser.parse_args()

    # Operations, outputs are read as an array of dicts to keep it ordered
    operations_array = []
    outputs_array = []

    try:

        for each in args.operation:

            # Split on the colon to specify the name on the left of it
            (name, instruction) = each.split(':')

            if instruction.find("<spacer>") > 0:
                print("Hey, you can't name a capture group "+
                    "'spacer', I'm using that!")
                exit(1)

            if instruction.find("<input>") > 0:
                print("Hey, you can't name a capture group "+
                    "'input', I'm using that!")
                exit(1)

            # similarly use the ` > ` to specify the flow of input to
            # the regex
            (input_string, regex_string) = re.split("\s>\s",
                instruction.strip())

            # append that to the operations array
            operations_array.append( [name, 
                    input_string.strip(), regex_string.strip()] 
                )

    except:

        # Failure likely from lack of operations to do
        print("\n"+"["+str(time.time())+"]"+" : "+"Wait a second, I don't understand the "+
            "operations to be done! Are there any? "+
            "Maybe there's small part I'm choking on? "+
            "Maybe try adding steps in one at a time in "+
            "interactive context with --limit set. "+
            "Exiting...")
        exit(1)

    try:

        for each_id, each_seq in zip(args.output_id, args.output_seq):

            # append that to the outputs array
            outputs_array.append( [each_id, each_seq] )

    except:

        # Failure likely from lack of operations to do
        print("\n"+"["+str(time.time())+"]"+" : "+
            "Wait a second, I don't understand the "+
            "outputs to be done! Are there any? "+
            "Maybe there's small part I'm choking on? "+
            "Maybe try adding steps in one at a time in "+
            "interactive context on a test set? "+
            "Exiting...")
        exit(1)

    if args.verbose > 0:
        print("\n"+"["+str(time.time())+"]"+" : "+"I'm reading in '"+#vars(args)["input-fastq"]+"', "+
            #"treating it as a "+("zipped" if args.gzipped else "unzipped")+" file and "+
            "applying these operations of alignment :\n")
        for each in operations_array:
            print("- "+each[0]+" :\n"+
                "  from : "+each[1]+"\n"+
                "  extract groups with regex : '"+each[2])

    if args.verbose > 0:
        print("\n"+"["+str(time.time())+"]"+" : "+"...and with these filters:")
        try:
            for i in args.filter:
                print("  "+i)
        except:
            print("  ( no filters defined )")
            args.filter = ["True == True"]

    if args.verbose > 0:
        print("\n"+"["+str(time.time())+"]"+" : "+"Then, I'm going to write out a FASTQ file to '"+
            "here or stdout")#vars(args)["output-base"]+".fastq'",end="")
        if args.report is not None:
            print(" and a report to '"+vars(args)["report"]+"'",end="")
        print(".")

#    # checking file existance for outputs, zipping together the 
#    # output base with each of the three. 
#    exit_flag = 0
#    for each in zip( [vars(args)["output-base"]]*20,            \
#                    ["_fail.fastq", "_pass.fastq",              \
#                        "_report.fastq", "_report.csv" ] ):
#        # At this stage, the tuple is joined to make the filename
#        this_path = ''.join(each)
#        # If the write-report flag is off and the path is the report,
#        # then this won't trip True for that path existing
#        if os.path.isfile(this_path) and              \
#                not( not(args.write_report) and       \
#                    (this_path.find("_report")>=0) ):
#            print("\n"+"["+str(time.time())+"]"+" : "+"File "+this_path+
#                " exits, so I'm quitting before you ask me to do "+
#                "something you might regret.")
#            exit_flag = 1
#    if exit_flag == 1:
#        exit(1)

#    # memory debugging
#    if args.memory_tracking:
#        tracemalloc.start(10)

    # We begin
    if args.verbose > 0:
        print("\n"+"["+str(time.time())+"]"+" : BEGIN\n")
    
    reader(
        vars(args)["input"],
        args.gzipped,
        vars(args)["output"],
        vars(args)["failed"],
        vars(args)["report"],
        args.verbose,
       # args.memory_tracking,
        operations_array, 
        args.filter ,
        outputs_array
        )

    print("\n"+"["+str(time.time())+"]"+" : "+
        "All worked 'till the work is done --- or some fatal error.")
    exit(0)
