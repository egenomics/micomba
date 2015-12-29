# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 16:43:43 2013

@author: jl
"""
from Bio import SeqIO
from functions import *

threads = 7 #This number should be the same as the number of cores used in the extractor.py
records = []
records_orig = []

#Gather the .fastq and .fa files generated by the different cores and put them together
i=1
while i<=threads:
    infile='output/output_'+str(i)+'.fa'
    infile2='output/output_'+str(i)+'.fastq'
    for record in SeqIO.parse(infile, "fasta"):
        records.append(record)
    for record in SeqIO.parse(infile2, "fastq"):
        records_orig.append(record)
    print i
    i+=1
SeqIO.write(records, "output/output_all.fasta", "fasta") #Export all the cut reads to a fasta file
SeqIO.write(records_orig, "output/output_all.fastq", "fastq") #Export all the original filtered reads to a  fastq file

#for each read, look at dictionary and redirect it to the proper file  
letters= map(chr, range(65, 73)) #=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

#Generate dictionary structure for individuals and locus
records={} #cut reads in fasta
records2={} #original filtered reads in fastq
for letter in letters:
    i = 1
    while i < 13:
        individual=letter+format(i)
        records[individual]={}
        records2[individual]={}
        i+=1

#Read output (cut fasta) and classify it into individual and locus dictionary structure
infile='output/output_all.fasta'
for record in SeqIO.parse(infile, "fasta"):
    this_record=record.id.split('_')
    if this_record[1] in records[this_record[0]]:
        records[this_record[0]][this_record[1]].append(record)
    else:
        records[this_record[0]][this_record[1]]=[]
        records[this_record[0]][this_record[1]].append(record)

#Read output (original filtered fastq) and classify it into individual and locus dictionary structure
infile='output/output_all.fastq'
for record in SeqIO.parse(infile, "fastq"):
    this_record=record.id.split('_')
    if this_record[1] in records2[this_record[0]]:
        records2[this_record[0]][this_record[1]].append(record)
    else:
        records2[this_record[0]][this_record[1]]=[]
        records2[this_record[0]][this_record[1]].append(record)

try:
    os.system('mkdir output')
except:
    pass

#We will parse the dictionary structure of records, navigating individuals and locus
for indiv in records:
    try:
        os.system('mkdir output/'+indiv) #Folder for individual
    except:
        pass
    for locus in records[indiv]:
        try:
            os.system('mkdir output/'+indiv+'/'+locus) #Folder for locus
        except:
            pass
        try:
            #Write output and run mira assembly
            SeqIO.write(records[indiv][locus], "output/"+indiv+'/'+locus+'/'+indiv+'_'+locus+'_output.fa', "fasta")
            SeqIO.write(records2[indiv][locus], "output/"+indiv+'/'+locus+'/'+indiv+'_'+locus+'_output.fastq', "fastq")
            if len(records[indiv][locus])>=5:
                #mira_assembly(indiv,locus)
                #cluster_and_align(indiv,locus,0.95)#individual, locus and clustering percentage
                f = open('output/'+indiv+'/'+locus+'/'+indiv+'_'+locus+'_lengths.txt','w')
                for record in records[indiv][locus]:
                    print>>f, record.id+'\t'+record.id.split('_')[3]
                f.close()
                R_plot(indiv,locus) #Generate R code
            else:
                print 'Not enough reads here', indiv, locus
        except:
            print 'No reads here', indiv, locus

