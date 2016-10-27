# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 14:08:14 2016

@author: karim
"""

#définir le directory:
import os 
os.chdir('C:\\Users\\karim\\FirstStepProject\\.spyproject')

#importer les séequence grace à l'URL des fichiers:
import requests
#la séquence complete
url = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21_GL383578v2_alt.fa.gz'
filename = url.split("/")[-1]
with open(filename, "wb") as f:
    r = requests.get(url)
    f.write(r.content)
    
#Décompresser la séquence et l'ouvrir
#décompresser le fichier  #la séquence est stoquée à présent dans la string file_content
import gzip
with gzip.open('C:\\Users\\karim\\FirstStepProject\\.spyproject\\chr21_GL383578v2_alt.fa.gz', 'rb') as f:
    file_content = f.read()
    
file_string = str(file_content,'utf-8') #type(file_string) Out[19]: str permet de passer d'un 8bytà un string simple (variante du meme objet mais plus manipulable)
file_string2 = file_string[file_string.find('\n')+1:] #la première ligne content l'id est supprimée
file_string3 = file_string2.replace('\n','') 
file_string4 = file_string3.lower()

#enregistrement de la séquences (string) dans un fichier texte pour simplifier la génération du pdf
text_file = open("Output.txt", "w")
text_file.write(file_string4)
text_file.close()

#insérer un retour à la ligne tous les n caractère --> adaptation à la page et la police utilisée par report lab
# Création d'un pdf avec numéro de pages
#conversion de la string en texte puis pdf
from reportlab.pdfgen import canvas
from reportlab.pdfgen.canvas import Canvas
c = canvas.Canvas('kemp.pdf')
y = 700
for line in open('Output.txt','r'):
    c.drawString(100, y, line.strip('\n'))
    c.showPage()
    c.save()
    