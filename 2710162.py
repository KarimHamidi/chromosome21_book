# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:52:11 2016

@author: karim
"""
#TEXT FILE OF THE CDS
#http://www.uniprot.org/docs/humchr21.txt

import os
os.chdir('C:\\Users\\karim\\FirstStepProject\\.spyproject')
import requests
list_url = ['http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21_GL383578v2_alt.fa.gz','http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21_GL383579v2_alt.fa.gz','http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21_GL383580v2_alt.fa.gz','http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21_GL383581v2_alt.fa.gz','http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21_KI270872v1_alt.fa.gz','http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21_KI270873v1_alt.fa.gz','http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21_KI270874v1_alt.fa.gz','http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz']

for i in list_url:
    filename = i.split("/")[-1]
    with open(filename, "wb") as f:
        r = requests.get(i)
        f.write(r.content)

#décompresser les fichiers et les sauver dans le meme directory   
import os      
import gzip
import glob
import os.path
source_dir = 'C:\\Users\\karim\\FirstStepProject\\.spyproject'
dest_dir = 'C:\\Users\\karim\\FirstStepProject\\.spyproject'

for src_name in glob.glob(os.path.join(source_dir, '*.gz')):
    base = os.path.basename(src_name)
    dest_name = os.path.join(dest_dir, base[:-3])
    with gzip.open(src_name, 'rb') as infile:
        with open(dest_name, 'wb') as outfile:
            for line in infile:
                outfile.write(line)

with open('chr21.fa', 'r') as infile:
    Séquence_string = ''
    for line in infile:
        if '>' not in line:
            Séquence_string = Séquence_string + str(line)
            Séquence21_string = Séquence_string.replace('N',' ').lower()    

with open('chr21_GL383578v2_alt.fa', 'r') as infile:
    Séquence_string = ''
    for line in infile:
        if '>' not in line:
            Séquence_string = Séquence_string + str(line)
            Séquence_Test_string = Séquence_string.replace('\n','').lower()              
               
#Génération d'un pdf Test
import textwrap
NewString = textwrap.fill(Séquence_Test_string, 110)

text_file = open("TEST.txt", "w")
text_file.write(NewString)
text_file.close()

#A présent la séquence est sous forme de Fichier Text et de string avec un retour à la ligne 
#tous les 50 caractères. 
#Génération des pages suivantes platypus? Text?
from reportlab.lib.enums import TA_JUSTIFY
from reportlab import canvas
title = 'Chromosme 21'
ptr = open("TEST.txt", "r")  #fichier text à convertir
lines = ptr.readlines()
ptr.close()
i = 800
numeroLine = 0
c = canvas.Canvas('Chr21C60999.pdf')
image = 'IMAGE.jpg'
c.setFont('Helvetica',60, leading= None)
c.drawCentredString(300, 800,title)
c.drawImage(image,100,100,width=200, height =200)
c.showPage()
while numeroLine < len(lines):
    if numeroLine - len(lines) < 60:
        # Défini les nombre de ligne par page
        i=800
        for line in lines[numeroLine:numeroLine+60]:
            c.setFont('Helvetica',10,leading=None)
            c.drawString(15, i, line.strip())   #si centré 300 sinon 15
            numeroLine += 1
            i -= 12                 #espace interligne
        c.showPage()
#    else:
#        i = 800
#        for line in lines[numeroLine:]:
#           numeroLine += 1
#           i -= 12

c.save()

#modifié celui d'au-dessus fonction ne pas toucher!!!!!
title = 'Chromosme 21'
title2 = 'Table des matières'
ptr = open("TEST.txt", "r")  #fichier text à convertir
lines = ptr.readlines()
ptr.close()
i = 800
longue = str(len(lines)*110) #à corriger
image = 'image2.png'
numeroLine = 0
c = canvas.Canvas('ChrPAGE.pdf')
#TITRE
c.setFont('Helvetica',60, leading= None)
c.drawCentredString(300, 650,title)
#Nombre de nucléotides annoncées.
c.setFont('Helvetica',30,leading=None)
c.drawCentredString(300,500,"Nombre totale de nucléotides")
c.drawCentredString(300,450,longue)
c.showPage()
#première image
c.setFont('Helvetica',45,leading=None)
c.drawCentredString(300,750,title2)
c.drawImage(image,80,200,width=300, height =400)
c.showPage()
#Séquence
while numeroLine < len(lines):
    if numeroLine - len(lines) < 60:
        # Défini les nombre de ligne par page
        i=800
        for line in lines[numeroLine:numeroLine+60]:
            c.setFont('Helvetica',10,leading=None)
            c.drawString(15, i, line.strip())   #si centré 300 sinon 15
            numeroLine += 1
            i -= 12             
            #espace interligne
    for i in range(10):
        c.setFont('Helvetica',8,leading=None)
        page_num = c.getPageNumber()
        text = "page %s" % page_num
        c.drawCentredString(300, 40, text)
    c.showPage()
#    else:
#        i = 800
#        for line in lines[numeroLine:]:
#           numeroLine += 1
#           i -= 12

c.save()



for i in range(10):
        page_num = c.getPageNumber()
        text = "This is page %s" % page_num
        c.drawString(100, 750, text)
        c.showPage()
    c.save()





