import sys
import os
import argparse
import time

from collections import defaultdict
import pylab as pl
import numpy as np
from Bio import SeqIO

def check_arguments():
    
    parser = argparse.ArgumentParser(usage='Script for analysis of fastq files')
    
    if(len(sys.argv) != 3):
          parser.error('''Invalid number of arguments!
                Needed arguments : <INPUT-FQ-FILE>  <OUTPUT-FOLDER>
                
                INPUT-FQ-FILE - fastq file that is going to be processed
                OUTPUT-FOLDER - directory in which the results are going to be stored''')
          exit()
    
    parser.add_argument("fastq", help="Fastq file that needs to be processed")
    parser.add_argument("outputDir", help="Directory in which the results are going to be stored.")
    
    args = parser.parse_args()

    if  not (args.fastq.endswith('.fastq') or args.fastq.endswith('.fq')) or (not os.path.isdir(args.outputDir)):
        parser.error('''Invalid types of given arguments!
                Needed arguments : <INPUT-FA-FILE>  <OUTPUT-FOLDER>
                
                INPUT-FQ-FILE - fastq file that is going to be processed
                OUTPUT-FOLDER - directory in which the results are going to be stored''')
        exit()
 
    
    
def makeHistogram(destinationFolder, lengthDict):
    pl.hist(lengthDict.keys(), len(lengthDict), weights=lengthDict.values())
    pl.title("Sequence length distribution")
    pl.savefig(destinationFolder + '/FastqLengthHist.png')



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
 
 
    
def statistics(fastqPath, destinationFile, numberOfReads, averageLength, averageReadQuality):
    fastqFile = open(fastqPath)
    statisticsFile = open(destinationFile + "/FastqStatistics.txt", "w+")
    statisticsFile.write("Fastq file : " + os.path.basename(fastqPath) + "\n")
    bytes = os.path.getsize(fastqPath)
    size = human_readable(bytes)
    statisticsFile.write("Fastq file size = " + str(size))
    statisticsFile.write("Number of reads = " + str(numberOfReads))
    statisticsFile.write("\nAverage read length = " + str(averageLength))
    statisticsFile.write("\nAverage read quality = " + str(averageReadQuality))
    statisticsFile.write("\nDate of analysis = " + str(time.strftime("%d/%m/%Y")))
    
    
def main():
    
    check_arguments()
     
    lenghtDict = defaultdict(int)
 
    numberOfReads = 0
    sumOfLengths = 0
    sumOfQuality = 0
    
    for record in SeqIO.parse(sys.argv[1], "fastq"):
        numberOfReads = numberOfReads + 1
        sumOfLengths = sumOfLengths + len(record)        
        lenghtDict[len(record.seq)] = lenghtDict[len(record.seq)] + 1
        sumOfQuality += sum(record.letter_annotations["phred_quality"]) / float(len(record.seq))
    
    makeHistogram(sys.argv[2], lenghtDict)
    statistics(sys.argv[1], sys.argv[2], numberOfReads, sumOfLengths / float(numberOfReads), sumOfQuality / numberOfReads)
 
if __name__ == '__main__':
    main()
    
