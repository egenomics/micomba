# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 11:13:08 2013

@author: JL Villanueva-Cañas
@email: jlvillanueva84@gmail.com

miraconvert -t hsnp A10_A04_first_assembly_out.maf assembly.html
"""

#Import all the functions and modules used. Check functions.py
import sys
from functions import *

try:
    os.system('mkdir output')
except:
    pass
#Open file containing tags (F* and R*)
records_filtered=[]
records_cut=[]
infiles=['IE7ZMMF01.fastq','IE7ZMMF02.fastq'] #455.094 total number of reads
thread = int(sys.argv[1]) #Current thread
threads = int(sys.argv[2]) #Number of threads/processors to use

#Asking for user input
tags = ['BarcodesNEW.fasta','Primers.fasta']
#infile = raw_input("Enter the path to the file containing "+bcolors.OKBLUE+"tags"+bcolors.ENDC+": ")
#infile2 = raw_input("Enter the path to the file containing "+bcolors.OKGREEN+"reads"+bcolors.ENDC+": ")
errors_tags=2
errors_primers=3

tag_list = get_tag_flavours(tags[0]) #Obtaining all tag variations
primers_list=get_primer_flavours(tags[1]) #Obtaining all primer variations
print "I'm looking for barcodes and primers...Get a coffee, this may take a while"

for infile in infiles:
    records = list(SeqIO.parse(infile, infile[infile.find('.')+1:])) #infile[infile.find('.')+1:] <--gets automatically the extension of the file
    start=len(records)/threads*(thread-1)
    if thread == threads:
        end=len(records)
    else:
        end=len(records)/threads*(thread)
    records = records[start:end]
    for record in records:
        barcode = find_pair(record, tag_list, errors_tags, thread) #Look for the barcodes
        if barcode != 'not_found':
            primers = find_pair(barcode, primers_list, errors_primers, thread)
            if primers!= 'not_found':
                record_formatted, record_cut=remake_id(primers, primers_list, thread)
                if record_formatted != 'locus_error':
                    records_filtered.append(record_formatted) #Save the read in a list
                    records_cut.append(record_cut) #Save the read in a list

SeqIO.write(records_filtered, "output/output"+'_'+format(thread)+'.fastq', "fastq") #Export all the selected reads to a file
SeqIO.write(records_cut, "output/output"+'_'+format(thread)+'.fa', "fasta") #Export all the selected reads to a file
