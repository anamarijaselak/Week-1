import sys
import os
import argparse
import math
import time

import pylab as pl
import numpy as np
from Bio import SeqIO
from _collections import defaultdict

def check_arguments():
    
    parser = argparse.ArgumentParser(usage='Script for analysis of fastq files')
    
    if(len(sys.argv) != 3):
          parser.error('''Invalid number of arguments!
                Needed arguments : <INPUT-FA-FILE>  <OUTPUT-FOLDER>
                
                INPUT-FA-FILE - fasta file that is going to be processed
                OUTPUT-FOLDER - directory in which the results are going to be stored''')
          exit()
    
    parser.add_argument("fasta", help="Fasta file that needs to be processed")
    parser.add_argument("outputDir", help="Directory in which the results are going to be stored.")
    
    args = parser.parse_args()

    if  not (args.fasta.endswith('.fasta') or args.fasta.endswith('.fa')) or (not os.path.isdir(args.outputDir)):
        parser.error('''Invalid types of given arguments!
                Needed arguments : <INPUT-FA-FILE>  <OUTPUT-FOLDER>
                
                INPUT-FA-FILE - fasta file that is going to be processed
                OUTPUT-FOLDER - directory in which the results are going to be stored''')
        exit()
       
def makeHistogram(destinationFolder, entropiesDict):
    pl.hist(entropiesDict.keys(), len(entropiesDict), weights=entropiesDict.values())
    pl.title("Sequence entropy distribution")
    pl.savefig(destinationFolder + '/FastaEntropiesHist.png')
    

def human_readable(bytes):
    size = ""
    if(bytes / pow(1000.0, 3) > 1):
        size = str(bytes / pow(1000, 3)) + " GB "
        bytes = bytes % pow(1000, 3)
    if(bytes / pow(1000.0, 2) > 1):
        size = size + str(bytes / pow(1000, 2)) + " MB "
        bytes = bytes % pow(1000, 2)
    if(bytes / 1000.0 > 1):
        size = size + str(bytes / 1000) + " KB "
        bytes = bytes % 1000
    if(bytes > 0):
        size = size + str(bytes) + " B "
    return size + "\n"
    
def statistics(fastaPath, destinationFile, numberOfReads, numberOfBadReads, numberOfKeptReads):
    fastaFile = open(fastaPath)
    statisticsFile = open(destinationFile + "/FastaStatistics.txt", "w+")
    statisticsFile.write("Fasta file : " + os.path.basename(fastaPath) + "\n")
    bytes = os.path.getsize(fastaPath)
    size = human_readable(bytes)
    statisticsFile.write("Fasta file size = " + str(size))
    statisticsFile.write("Number of reads = " + str(numberOfReads))
    statisticsFile.write("\nNumber of kept reads= " + str(numberOfKeptReads))
    statisticsFile.write("\nNumber of bad reads = " + str(numberOfBadReads))
    statisticsFile.write("\nDate of analysis = " + str(time.strftime("%d/%m/%Y")))
 
    
def main():

    check_arguments()
    
    entropiesDict = defaultdict(int)
    max_entropy = 0.
    numberOfReads = 0
    numberOfKeptReads = 0
    numberOfBadReads = 0
    newFastaFile = open(sys.argv[2] + "/newFastafile.fa", "w+")
    for record in SeqIO.parse(sys.argv[1], "fasta"):
        numberOfReads = numberOfReads + 1
        basesInRecord = 0
        baseDict = defaultdict(int)
        if "N" in record.seq :
            numberOfBadReads = numberOfBadReads + 1
            continue
        entropy = 0
        
        for base in record.seq:
            baseDict[base] = baseDict[base] + 1
            basesInRecord = basesInRecord + 1
            
        for key in baseDict.keys():
            entropy = entropy + (baseDict[key] / float(basesInRecord)) * math.log(baseDict[key] / float(basesInRecord), 2)
        
        entropiesDict[round(-1 * entropy, 2)] = entropiesDict[round(-1 * entropy, 2)] + 1
    
        if -1 * entropy > 0.5:
            newFastaFile.write(">")
            newFastaFile.write(record.name + "\n")
            newFastaFile.write(str(record.seq) + "\n")
            numberOfKeptReads = numberOfKeptReads + 1
        else :
            numberOfBadReads = numberOfBadReads + 1
    
    makeHistogram(sys.argv[2], entropiesDict)
    statistics(sys.argv[1], sys.argv[2], numberOfReads, numberOfBadReads, numberOfKeptReads)
 
   
        
if __name__ == '__main__':
    main() 
        
    
    
    
    
    
