# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 15:29:38 2016

@author: karim
"""
import os 
os.chdir('C:\\Users\\karim\\FirstStepProject\\.spyproject')

import requests

url = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21_GL383578v2_alt.fa.gz'
filename = url.split("/")[-1]
with open(filename, "wb") as f:
    r = requests.get(url)
    f.write(r.content)
    
import gzip

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio import SeqIO
for seq_record in SeqIO.parse("chr21_GL383578v2_alt.fa", "fasta"):
    print(seq_record.id) #nom de la séquence
    print(repr(seq_record.seq)) #séquence sans le nom, juste les nucléotides
    print(len(seq_record)) #longueur de la séquence

print(seq_record.seq)
seq1 = seq_record.seq
import reportlab
from reportlab.pdfgen import canvas
def hello(c):
    c.drawString(100,100,"essai")
c = canvas.Canvas("hell.pdf")
hello(c)
c.showPage()
c.save()

#création d'un fichier et compilation de toutes les séquences.
open("compilationtotal", 'a')  #création du fichier compilation
filenames = ['chr21_GL383578v2_alt.fa', 'chr21_GL383579v2_alt.fa', 'chr21_GL383580v2_alt.fa', 'chr21_GL383581v2_alt.fa', 'chr21_KI270872v1_alt.fa', 'chr21_KI270873v1_alt.fa', 'chr21_KI270874v1_alt.fa ']
with open('C:\\Users\\karim\\FirstStepProject\\.spyproject\\compilationtotal', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

#maintenant je cherche à éliminer les lignes contenant les id des contigs. Elles commencent par ">"

#for i, line in enumerate("compilationtotale"):
#        if not line.startswith('>'):
#            output.write(line)
#            
#for line in compilationtotale:
#    if not line.startswith('>'):
#        with open('C:\\Users\\karim\\FirstStepProject\\.spyproject\\compilationtotal', 'w') as outfile:
            
#AMELIORATION MTN LES FICHIER DES CONTIGS SONT COMPILES SANS LES ID DES SEQUENCES            
open("seqtotal", 'a')  #création du fichier compilation
filenames = ['chr21_GL383578v2_alt.fa', 'chr21_GL383579v2_alt.fa', 'chr21_GL383580v2_alt.fa', 'chr21_GL383581v2_alt.fa', 'chr21_KI270872v1_alt.fa', 'chr21_KI270873v1_alt.fa', 'chr21_KI270874v1_alt.fa ']
with open('C:\\Users\\karim\\FirstStepProject\\.spyproject\\seqtotal', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                if not line.startswith('>'):
                    outfile.write(line)           
            
#MAINTENANT IL FAUT COMPARER LES FICHIERS --> ESSAI QUI CONTIENT LES 7 CONTIGS ET LE FICHIER DE LA SEQ COMPLETE
#mais d'abord essai de changer les majuscule en minuscule
from itertools import chain
from glob import glob

file = open('seqtotal', 'r')

lines = [line.lower() for line in file]
with open('seqtotal', 'w') as out:
     out.writelines(lines)

#Traiter le fichier de la sequence complete de la meme manère!  
#1 supprimer l'id de la séquence 
open("chr22", 'a')

filenames = ['chr212014.fa ']
with open('C:\\Users\\karim\\FirstStepProject\\.spyproject\\chr22', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                if not line.startswith('>'):
                    outfile.write(line)    
#tout mettre en minuscule pour le fichier de la sequence complète du chr21
from itertools import chain
from glob import glob

file = open('chr22', 'r')

lines = [line.lower() for line in file]
with open('chr22', 'w') as out:
     out.writelines(lines)
#MAINTENANT LE FICHIER CONTENANT LA TOTALITE DE LA SEQUENCE NE CONTIENT PLUS QUE DES MINUSCULE LUI AUSSI

#compter les n et toutes les nucleotides par la meme occasionc
#import string 
#text = open('chr22')
#letters = string.ascii_lowercase
#for i in text:
#  text_lower = i.lower()
#  text_nospace = text_lower.replace(" ", "")
#  text_nopunctuation = text_nospace.strip(string.punctuation)
#  for a in letters:
#    if a in text_nopunctuation:
#      num = text_nopunctuation.count(a)
#      print(a, num) 
import gzip
with gzip.open('/home/joe/file.txt.gz', 'rb') as f:
    file_content = f.read()

with open('chr21_GL383578v2_alt.fa', 'r') as myfile:
    data=myfile.read().replace('\n', '')
    
print(len(data))

data.count('A')