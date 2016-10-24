# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:57:41 2016

@author: karim
"""
import os 
os.chdir('C:\\Users\\karim\\FirstStepProject\\.spyproject')
#from reportlab.platypus import image
import requests
import gzip
url = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21_GL383578v2_alt.fa.gz'
filename = url.split("/")[-1]
with open(filename, "wb") as f:
    r = requests.get(url)
    f.write(r.content)
    
#Décompresser la séquence et l'ouvrir
#décompresser le fichier  #la séquence est stoquée à présent dans la string file_content
with gzip.open('C:\\Users\\karim\\FirstStepProject\\.spyproject\\chr21_GL383578v2_alt.fa.gz', 'rb') as f:
    file_content = f.read()
    
#passer par une string pour créer un fichier texte. Le fichier sequence.txt contient la séquence    
import textwrap 
file_string = str(file_content,'utf-8') #type(file_string) Out[19]: str permet de passer d'un 8bytà un string simple (variante du meme objet mais plus manipulable)
file_string2 = file_string[file_string.find('\n')+1:] #la première ligne content l'id est supprimée
file_string3 = file_string2.lower()
file_string4 = file_string3.replace('\n','')
file_string5 = textwrap.fill(file_string4,80)

#insérer la séquence dans un fichier texte                           
text_file = open("Seqmin80.txt", "w")
text_file.write(file_string5)
text_file.close()

#GENERATION DU PDF :
from reportlab import *
#Génération des pages suivantes platypus? Text?
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
from PIL import Image
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
 

#page de titre basique. 
def generate_(c):
    image = 'Human_chromosome_21_from_Gene_Gateway_-_no_label.png'
    c.setFont('Helvetica',48, leading = None)
    c.drawCentredString(350,600,'Chromosome 21')
    c.drawImage(image, 100, 100,width = 200, height = 200)
    c.showPage()
#commande d'exécution de la fonction       
c = canvas.Canvas('21.pdf', pagesize = letter)
generate_book(c)    

#génération d'un pdf contenant la séquence

ptr = open("Seqmin80.txt", "r")  #fichier text à convertir
lines = ptr.readlines()
ptr.close()
i = 800
numeroLine = 0
c = canvas.Canvas('Chr21C.pdf')


while numeroLine < len(lines):
    if numeroLine - len(lines) < 63: # Défini les nombre de ligne par page
        i=800
        for line in lines[numeroLine:numeroLine+63]:
            c.setFont('Helvetica',10,leading=None)
            c.drawString(15, i, line.strip())   #si centré 300 sinon 15
            numeroLine += 1
            i -= 12                 #espace interligne
        c.showPage()
    else:
        i = 800
        for line in lines[numeroLine:]:
           numeroLine += 1
           i -= 12
c.showPage()
c.save()



    