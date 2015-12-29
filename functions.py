# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 12:44:56 2013

@author: JL Villanueva-Cañas
@email: jlvillanueva84@gmail.com
"""

import re, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from fuzzysearch import * # We are implementing the Levenshtein distance here. 

def get_tag_flavours(infile): #From a set of tags it creates all possible variations and stores it in a file and in a list object
    handle = open(infile, "rU")
    tag_list = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    try:
        os.system('rm '+infile+'_flavours')
    except:
        pass
    f = open(infile+'_flavours', 'a')
    for tag in tag_list:
        print>>f, '>'+tag_list[tag].id     
        print>>f, tag_list[tag].seq
        print>>f, '>'+tag_list[tag].id+'C'
        print>>f, tag_list[tag].seq.complement()
        print>>f, '>'+tag_list[tag].id+'R'
        print>>f, tag_list[tag].seq[::-1]
        print>>f, '>'+tag_list[tag].id+'RC'
        print>>f, tag_list[tag].seq.reverse_complement()
    f.close()
    handle = open(infile+'_flavours', "rU")
    tag_list = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return tag_list

def get_primer_flavours(infile): #From a set of tags it creates all possible variations
    handle = open(infile, "rU")
    tag_list = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    try:
        os.system('rm '+infile+'_flavours')
    except:
        pass
    f = open(infile+'_flavours', 'a')
    for tag in tag_list:
        print>>f, '>'+tag_list[tag].id     
        print>>f, tag_list[tag].seq
        print>>f, '>'+tag_list[tag].id[:tag_list[tag].id.find('_')]+'C'+'_'+tag_list[tag].id[tag_list[tag].id.find('_')+1:]
        print>>f, tag_list[tag].seq.complement()
        print>>f, '>'+tag_list[tag].id[:tag_list[tag].id.find('_')]+'R'+'_'+tag_list[tag].id[tag_list[tag].id.find('_')+1:]
        print>>f, tag_list[tag].seq[::-1]
        print>>f, '>'+tag_list[tag].id[:tag_list[tag].id.find('_')]+'RC'+'_'+tag_list[tag].id[tag_list[tag].id.find('_')+1:]
        print>>f, tag_list[tag].seq.reverse_complement()
    f.close()
    handle = open(infile+'_flavours', "rU")
    tag_list = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return tag_list

class bcolors: #Silly function to output colors and bold text
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''
        
def find_tag(pattern,text,errors):
    if errors==0:    
        match = re.search(pattern,text)
        if match: #First we look for exact match
            return [[match.start(),0,0]]
        else:
            return 'not_found'
    elif errors>0:
        match = find_near_matches_with_ngrams(pattern, text, errors)
        return match
            
def find_pair(record, tag_list, errors, thread):
    Forward=False #Forward and Reverse variables. They will turn True when we find a tag.
    Reverse=False
    f_mismatch=0 #Variables (forward and reverse) storing the Levenshtein distance. AKA the number of errors/differences between the tag and the read
    r_mismatch=0
    f_pos=0 #Position of forward tag
    r_pos=0 #Position of reverse tag
    error=0
    while error<=errors: #We'll start looking for 0 errors        
        for tag in tag_list:
            if tag[0]=='F' and Forward==False: #If we still don't have a hit for Forward
                match=find_tag(str(tag_list[tag].seq), str(record.seq), error)
                if match!='not_found': #We have a hit! Storing information
                    f_mismatch=match[0][2]
                    f_pos=match[0][0]
                    Forward=True
                    barcode1=tag[0:tag.find('_')]
            if tag[0]=='R' and Reverse==False: #Same for Reverse tag
                match=find_tag(str(tag_list[tag].seq), str(record.seq), error)
                if match!='not_found': #We have a hit! Storing information
                    r_mismatch=match[0][2]
                    r_pos=match[0][0]
                    Reverse=True
                    barcode2=tag
        error+=1
    if Forward==True and Reverse==True: #We found both forward and reverse tag
        if len(record.id.split('_'))==1:#ID remains unmodified, so we're dealing with barcodes and not primers yet        
            f = open('output/barcodes_used'+'_'+str(thread)+'.txt', 'a')
            record.id=record.id+'_'+barcode1+'_'+barcode2+'_'+format(f_mismatch)+'_'+format(r_mismatch) #Store in the ID the tags found
            print>>f, ('\t'.join(map(format,record.id.rstrip().split('_'))))
            f.close()
        else: #We are dealing with primers
            record.id=record.id+'_'+barcode1+'_'+barcode2+'_'+format(f_mismatch)+'_'+format(r_mismatch)+'_'+format(f_pos)+'_'+format(r_pos)+'_'+format(len(record)) #Store in the ID the tags found
        return record
    else: #generating output for the different scenarios where we don't find all elements
        if Forward==True:
            f = open('output/reads_one_barcode_primer'+'_'+str(thread)+'.txt', 'a')
            print>>f, record.id, barcode1, f_mismatch, f_pos
            f.close()
        elif Reverse==True:
            f = open('output/reads_one_barcode_primer'+'_'+str(thread)+'.txt', 'a')
            print>>f, record.id, barcode2, r_mismatch, r_pos
            f.close()
        if Forward==False and Reverse==False and len(record.id.split('_'))==1:
            f = open('output/reads_no_barcode'+'_'+str(thread)+'.txt', 'a')
            print>>f, record.id
            f.close()
        elif Forward==False and Reverse==False and len(record.id.split('_'))!=1:
            f = open('output/reads_barcode_no_primers'+'_'+str(thread)+'.txt', 'a')
            print>>f, record.id
            f.close()
        return 'not_found'

def remake_id(primers, primers_list, thread):
    cut_primers = SeqRecord(Seq(""),id="")
    f = open ( 'ID_tags.csv' , 'r')
    l = [ map(str,line.rstrip().split('\t')) for line in f ]
    tags_dict ={} #Creating dictionary from csv matrix    
    for item in l:
        i=1
        while i<len(l[0]):
            tags_dict[item[i]]=l[0][i]+item[0]
            i+=1
    tag_row = primers.id.split('_')
    tag_row.append(tag_row[1]+'_'+tag_row[2]) #Concatenate forward and reverse barcode
    nums = map(int, re.findall(r'\d+', tag_row[-1]))
    tag_row.append(tags_dict['F'+format(nums[0])+'_R'+format(nums[1])]) #Add individual information using barcodes
    if int(tag_row[-4])>=int(tag_row[-5]): #Forward primer first
        if tag_row[5]+'_'+tag_row[7] in primers_list:
            distance=int(tag_row[-4])-int(tag_row[-5])+len(primers_list[tag_row[6]+'_'+tag_row[7]].seq)
            cut_primers.seq=primers.seq[int(tag_row[-5]):int(tag_row[-4])+len(primers_list[tag_row[6]+'_'+tag_row[7]].seq)] #Select only sequence between primers
        else:
            f = open('output/primers_locus_no_match'+'_'+str(thread)+'.txt', 'a')
            print>>f, ('\t'.join(map(format,tag_row)))
            f.close()
            return ('locus_error','locus_error')
    else:
        if tag_row[5]+'_'+tag_row[7] in primers_list: #Reverse primer first
            distance=int(tag_row[-5])-int(tag_row[-4])+len(primers_list[tag_row[5]+'_'+tag_row[7]].seq)
            cut_primers.seq=primers.seq[int(tag_row[-4]):int(tag_row[-5])+len(primers_list[tag_row[5]+'_'+tag_row[7]].seq)] #Select only sequence between primers
            #print distance
        else:
            f = open('output/primers_locus_no_match'+'_'+str(thread)+'.txt', 'a')
            print>>f, ('\t'.join(map(format,tag_row)))
            f.close()
            return ('locus_error','locus_error')
    if 'RC' in tag_row[5]: #Reversing all the reads that contain an RC forward primer
        cut_primers.seq = cut_primers.seq[::-1].complement()
    tag_row.append(distance)
    primers.id = tag_row[14]+'_'+tag_row[7]+'_'+tag_row[0]+'_'+format(distance) #Individual_locus_read
    cut_primers.id = primers.id
    f = open('output/barcodes_primers_used'+'_'+str(thread)+'.txt', 'a')
    print>>f, ('\t'.join(map(format,tag_row)))
    f.close()
    return (primers, cut_primers)

def mira_assembly(indiv,locus,cluster_perc):
    original_dir = os.getcwd()
    os.chdir("output/"+indiv+'/'+locus+'/')
    f = open('manifest.conf', 'w')
    print>>f,'project = '+indiv+'_'+locus+'_first_assembly'
    print>>f,'job = genome,denovo,accurate'
    print>>f,'readgroup = Ind'+indiv+'_locus_'+locus
    print>>f,'data = '+indiv+'_'+locus+'_output.fastq'
    print>>f,'technology = 454'
    print>>f,'parameters = 454_SETTINGS -ASSEMBLY:minimum_reads_per_contig ≥ 2'
    f.close()
    os.system('mira manifest.conf > log_assembly.txt')
    os.chdir(original_dir)


def cluster_and_align(indiv,locus):
    il_dir="output/"+indiv+'/'+locus+'/'
    input_file = il_dir + +indiv+'_'+locus+'_output.fastq'
    output_file = il_dir + +indiv+'_'+locus+'_clustered.fastq'
    cmd='./cd-hit-v4.5.4-2011-03-07/cd-hit-454 -i ' + input_file + ' -o ' + output_file + ' -c  '+cluster_perc+  '-g  1 -D 1000';#cluster sequences together in target file  
    os.system(cmd)




def R_plot(indiv, locus):
    f = open('output/'+indiv+'/'+locus+'/'+indiv+'_'+locus+'_lengths.R','w')
    print>>f,'library(ggplot2)'
    print>>f,"df <- read.delim(file='output/"+indiv+'/'+locus+'/'+indiv+'_'+locus+"_lengths.txt', header=F) #First column"
    print>>f,"names(df)<-c('read','length')"    
    print>>f,'or<-min(df$length)-0.5'
    print>>f,"png('output/"+indiv+'/'+locus+'/'+indiv+'_'+locus+"_lengths.png', width=7+2/3, height=6+2/3, units='in', res=600)"
    print>>f,"ggplot(df, aes(x=length)) + geom_histogram(colour='black', fill='turquoise3', binwidth = 1, origin=or) + scale_x_continuous(breaks=round(seq(min(df$length), max(df$length), by =1),1)) + labs(x='Primer-to-primer sequence length', y='Number of reads', title='Locus "+locus+", individual "+indiv+"') + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))"
    print>>f,"dev.off()"
    f.close()
    os.system('R CMD BATCH output/'+indiv+'/'+locus+'/'+indiv+'_'+locus+'_lengths.R')
    
def generate_html(records, all_locus):
    f = open('output/overview_table.html','w')
    print>>f,'<HTML>'
    print>>f,'<script type="text/javascript" charset="utf-8" src="../DataTables/media/js/jquery.js"></script>'
    print>>f,'<script type="text/javascript" charset="utf-8" src="../DataTables/media/js/jquery.dataTables.js"></script>'
    print>>f,'<script type="text/javascript" charset="utf-8" src="../DataTables/media/js/num-html.js"></script>'
    print>>f,'<script type="text/javascript" charset="utf-8" src="../DataTables/media/js/FixedColumns.js"></script>'
    print>>f,'<style type="text/css" title="currentStyle">'
    print>>f,'    @import "../DataTables/media/css/demo_table.css";'
    print>>f,'    @import "../DataTables/media/css/demo_page.css";'
    print>>f,'</style>'
    print>>f,'<style>'
    print>>f,'a:link {color:#0099FF;text-decoration: none}    /* unvisited link */'
    print>>f,'a:visited {color:#CC66CC;text-decoration: none} /* visited link */'
    print>>f,'a:hover {color:#FF3333;text-decoration: none}   /* mouse over link */'
    print>>f,'a:active {color:#009900;text-decoration: none}  /* selected link */'
    print>>f,'</style>'
    print>>f,'<script type="text/javascript" charset="utf-8">'
    print>>f,'$(document).ready( function () {'
    print>>f,"    var oTable = $('#overview').dataTable( {"
    print>>f,'        "sScrollY": "800px",'
    print>>f,'        "sScrollX": "100%",'
    print>>f,'        "sScrollXInner": "100%",'
    print>>f,'        "bScrollCollapse": true,'
    print>>f,'        "bPaginate": false'
    print>>f,'    } );'
    print>>f,'    new FixedColumns( oTable );'
    print>>f,'} );'
    print>>f,'  </script>'
    print>>f,'<table id="overview" class="display">'   
    print>>f,'  <thead>'
    print>>f,'      <tr>'
    print>>f,'      <th>Individual</th>'
    for locus in all_locus:
        print>>f,'                <th>'+locus+'</th>'
    print>>f,'      </tr>'
    print>>f,'  </thead>'
    print>>f,'  <tbody>'
    for indiv in records:
        print>>f,'          <tr class="gradeC">'
        print>>f,'              <td>'+indiv+'</td>'
        for locus in all_locus:
            if locus in records[indiv]:
                if len(records[indiv][locus])>4:
                    print>>f,"      <td><A HREF='"+indiv+'/'+locus+'/'+indiv+'_'+locus+"_lengths.png'>"+format(len(records[indiv][locus]))+'</A></td>'
                else:
                    print>>f,'      <td>'+format(len(records[indiv][locus]))+'</td>'
            else:
                print>>f,'      <td>0</td>'
        print>>f,'          </tr>'
    print>>f,'            </tbody>'
    print>>f,'            </table>'
    print>>f,'            </HTML>'
    f.close()