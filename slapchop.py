#!/usr/bin/env python3

# Simply Looking At Pairwise Comparisons to Optimally Parition
#
# Modifying to be more general.
#

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

import re
import multiprocessing 
import itertools
import argparse
import subprocess
import sys
import time

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
      #  print(pairwise2.format_alignment(*aln[0]))
      #  print(aln[0])
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
      #  print(pairwise2.format_alignment(*aln[0]))
      #  print(aln[0])
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
      #  print(pairwise2.format_alignment(*aln[0]))
      #  print(aln[0])
        ScoreUMI = round(alnUMI[0][2],3)
        AlignmentStartUMI = alnUMI[0][3]
        AlignmentEndUMI = alnUMI[0][4]
    
        umiSeq = tailSeq[(alnUMI[0][3]):(alnUMI[0][4])]
      #  umiSeq = umiSeq[16:31:1]
      #  umiSeq = umiSeq[17:31:2]
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
#    inputFastqList = f.readlines()
    iterf = iter(f)
    while True:
      i = 0
      thisBite = []
      while i < biteSize:
#        anEntry = inputFastqList[0:4]
#        del inputFastqList[0:4]
#        if anEntry:
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

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=""+
    "slapChopper.py")
  parser.add_argument("--inputFastq",required=True)
  parser.add_argument("--processes",default=1)
  parser.add_argument("--biteSize",
    help="The size of bites to chomp off of the input file for multi-process",
    default=1000)
  parser.add_argument("--fixed1pattern")
  parser.add_argument("--fixed2pattern")
  parser.add_argument("--umipattern")
  parser.add_argument("--matchScore",default=2)
  parser.add_argument("--mismatchScore",default=0.01)
  parser.add_argument("--readGapOpenScore",default=-1)
  parser.add_argument("--readGapExtendScore",default=-1)
  parser.add_argument("--templateGapOpenScore",default=-1)
  parser.add_argument("--templateGapExtendScore",default=-1)
  parser.add_argument("--filters")
  parser.add_argument("--outputBase")
  parser.add_argument("--outputReport")
  parser.add_argument("--logFile")
  parser.add_argument("--maxQueueSize",help="in gigs",default=10)
  args=parser.parse_args()
  print("===")
  
  # FIRST, where are we writing to?
  print("Looking at "+args.inputFastq+", aligning these to "+
    args.fixed1pattern+" , "+
    args.fixed2pattern+" , "+
    args.umipattern+" , "+
    "with parameters "+
    "Match "+str(args.matchScore)+", Mismatch "+str(args.mismatchScore)+
    ", ReadGapOpen "+str(args.readGapOpenScore)+
    ", ReadGapExtend "+str(args.readGapExtendScore)+
    ", TemplateGapOpen "+str(args.templateGapOpenScore)+
    ", TemplateGapExtend "+str(args.templateGapExtendScore)+
    "\n"
  )
  
  try:
    filterz = re.split(",",args.filters)
    print("Applying filters of "+str(filterz))
  except:
    None
  try:
    print("Writing report to "+args.outputReport)
  except:
    None
  try:
    print("Writing output to "+args.outputBase,end='')
    if args.filters:
      print("_[(pass)|(fail)].fastq")
    else:
      print("_all.fastq\n")
  except:
    None
  try:
    print("Doing this in bites of "+args.biteSize)
  except:
    None
  try:
    print("Writing progress to "+args.logFile)
  except:
    None
  
  print("BEGIN")
  print()
  
#  subprocess.run(" rm -f "+args.outputBase,shell=True)
  subprocess.run(" rm -f "+args.outputBase+"_all.fastq",shell=True)
  subprocess.run(" rm -f "+args.outputBase+"_fail.fastq",shell=True)
  subprocess.run(" rm -f "+args.outputBase+"_pass.fastq",shell=True)
  subprocess.run(" rm -f "+args.outputReport,shell=True)
  try:
    subprocess.run(" rm -f "+args.logFile,shell=True)
  except:
    print("You're not using a log file???")

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
