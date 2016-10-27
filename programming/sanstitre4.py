## -*- coding: 3-8 -*-
#"""
#Created on Mon Oct  3 10:14:48 2016
#
#@author: karim
#"""
##EXERCICES 1  031016 6Q
#fasta=""">gi|506953611|gb|KC684925.1| Amphiprion clarkii rhodopsin (RH) mRNA, partial cds
#AGTCCTTATGAGTACCCTCAGTACTACCTTGTCAACCCAGCCGCTTATGCTGCTCTGGGTGCCTACATGT
#TCTTCCTCATCCTTGCTGGCTTCCCAGTCAACTTCCTCACCCTCTACGTCACCCTCGAACACAAGAAGCT
#GCGAACCCCTCTAAACTACATCCTGCTGAACCTCGCGGTGGCTAACCTCTTCATGGTGCTTGGAGGATTC
#ACCACAACGATGTACACCTCTATGCACGGCTACTTCGTCCTTGGACGCCTCGGCTGCAATCTGGAAGGAT
#TCTTTGCTACCCTCGGTGGTGAGATTGCCCTCTGGTCACTGGTTGTTCTGGCTATTGAAAGGTGGGTCGT
#TGTCTGCAAGCCCATCAGCAACTTCCGCTTCGGGGAGAATCACGCTATTATGGGTTTGGCCTTCACCTGG
#ACAATGGCCAGTGCCTGCGCTGTTCCTCCTCTTGTCGGCTGGTCTCGTTACATCCCTGAGGGCATGCAGT
#GCTCATGTGGAGTTGACTACTACACACGTGCAGAGGGTTTCAACAATGAGAGCTTTGTCGTCTCCTCTTG
#TCGGCTGGTCTCGTTACATCCCTGAGG"""
#
##Q1
##le plus grand que indique le début de la séquence. Il faut pouvoir lire n'importe quelle fichier fasta.
##la première ligne est constante dans un fichier fasta. donc il faut trouver comment compter juste la seq. 
##ici on ne supprime pas la première ligne d'identité
##la séquence commence à la ligne suivante. Attention la séquence peut etre précédée de tabulateurs ou de retour
##la première étape: commencer par la longueur totale? c'est une possibilité
##puis soustraire la première ligne /n (position du premier /n) ensuite il faut compter les slash n et soustraire encore
##sinon #la commande ne passe pas
#start=fasta.find('\n')+1
#seq=fasta[start:].replace(" ","").replace("\n","")
#lenDNA=len(seq)
#print("length of the sequence =", lenDNA)
#
##question 2 calculer la fréquence des bases
#
#print(seq)
#fA=seq.count("A")/lenDNA
#fC=seq.count("C")/lenDNA
#fG=seq.count("G")/lenDNA
#fT=seq.count("T")/lenDNA
#print("frequency of A, C, G, T, = %f, %f, %f, %f" % (fA,fC,fG,fT))
##Question 3
#first=seq[0::3]
#second=seq[1::3]
#third=seq[2::3]
#print("First codon position", first)
#print("Second codon position", second)
#print("third codon position", third)
##Question 4
#gbnum=fasta.split('|')[3]
#print("Genbak Nubmer:", gbnum)
##Question 5 
#print("mRNA sequece:", seq.replace("T","U"))
##Question 6 Combien d'acides aminés
#print("Number of Amino Acids:", lenDNA//3)
#fin de l'exercice 1
#Tupple
#EXERCICES SUR LES DICTIONNAIRES. 
#Question 1
#for seq in (seq1, seq2,)

#dico={}
##for seq in (seq1, seq2)
#seq1=">Amphiprion clarkii\nAGTTGACCTAGTCATAGA"
#seq2=">Amphiprion frenatus\nAGCTGACCTAGTTTTAGA"
#seq3=">Amphiprion ocellaris\nAGTTGACCTGGGCATCGA"
#seq4=">Pomacentrus mollucensis\nAGTCTACCTGATCCGGA"
#liste =seq1.split('\n')
#liste[0]=liste[0].strip('>')
#dico[liste[0]]=liste [1]
#
#liste =seq2.split('\n')
#liste[0]=liste[0].strip('>')
#dico[liste[0]]=liste [1]
#
#liste =seq3.split('\n')
#liste[0]=liste[0].strip('>')
#dico[liste[0]]=liste [1]
#
#liste =seq4.split('\n')
#liste[0]=liste[0].strip('>')
#dico[liste[0]]=liste [1]
#
#print(dico)
#
#seq1b="ATAATATTCGATTGATCAGT"
#seq2b="ATAATACTCGATTTATCAGT"
#seq3b="ATAATACTCGATCGATCCGT"
#seq4b="ATAATAGGCGATCGACTAGT"
#
#liste =seq1b.split('\n')
#liste[0]=liste[0].strip('>')
##keep the keys and replace the values
#dico2= dico.copy()
#dico2["Amphiprion clarkii"]=seq1b
#dico2["Amphiprion frenatus"]=seq2b
#dico2["Amphiprion ocellaris"]=seq3b
#dico2["Pomacentrus mollucensis"]=seq4b
#
#seq1g = seq1.count("C")
#seq1c = seq1.count("G")
#seq1gc= seq1g + seq1c 
#print("seq1 G+C:", seq1gc/ len(seq1))
#
#seq2g = seq2.count("C")
#seq2c = seq2.count("G")
#seq2gc= seq2g + seq2c
#print("seq2 G+C:",seq2gc/ len(seq2))
#
#seq3g = seq3.count("C")
#seq3c = seq3.count("G")
#seq3gc= seq3g + seq3c
#print("seq3 G+C:",seq3gc/ len(seq3))
#
#seq4g = seq4.count("C")
#seq4c = seq1.count("G")
#seq4gc= seq4g + seq4c
#print("seq4 G+C:",seq4gc/ len(seq4))
#
#GC1 = (seq1gc/ len(seq1))
#GC2 = (seq2gc/ len(seq2))
#GC3 = (seq3gc/ len(seq3))
#GC4 = (seq4gc/ len(seq4))
#
#GC =[GC1,GC2,GC3,GC4]
#GC.sort()
#
#k=list(dico.keys())
#
#gcitems = [(GC[0],k[0]),(GC[1],k[1]),(GC[2],k[2]),(GC[3],k[3])
#
#gcitems.sort() 

#z = 1 and 0
#print(z)

#EXERCICES 10.10.16
#ATTENTION REGARDER LA LONGUEUR DE LA SEQUENCE
#haveIt=False 
#x=input("Enter an Amino Acid") #OK
#seq="NFYIPMFNKTGVVRSPFEYPQYYLAGVVRSPFEY"
#for i in seq: 
#    if x == i:
#        haveIt = True 
#        
#        break
#if haveIt:
#    print(x, " is  in the sequnence")
#else:
#    print(x,"is not in the sequence")

#Question 2
#import random 
#dna = ["A","G","C","T"]
##initialise empty string
##this is where the bases will be added as they are generated
#random_sequence=''
## We'll create a string of 150 random bases
#for i in range(0,150):
#    random_sequence+=random.choice(dna)
#print(random_sequence)
#ATTENTION AUX ETAPES VOIR PDF
#QUESTION * 

#seq="NFYIPMFNKTGVVRSPFEYPQYYLAGVVRSPFEY"
#seq2="ALILALSMGY"
#count= 0
#for i in range(len(seq2)):
#    already=False
#    for k in range(i):
#        if seq2[k] == seq2[i]:
#            already=True
#    if not already and seq2[i] in seq:
#        count+=1
#print("there are %d aa from seq2 in seq" %count)
#
##plus simple c'est d'utiliser la fonction 
#intersect = set(seq2) & set(seq)
#print("there are %d AA from seq2 in seq"%len(intersect))



#24.10.2016 Exercices traitement de données
import requests
url = 'http://www2.unil.ch/phylo/teaching/python/clownfish.gb'
filename = url.split("/")[-1]
with open(filename, "wb") as f:
    r = requests.get(url)
    f.write(r.content)

#créer un dictionnaire vide
import string
dico ={}
noms = []
séquence=[]
f = open('clownfish.fasta','+r')
espece = 0
lines = f.readline()

if line != '/n':
    for line in f:
        if '>' in line:
            noms.append(line)
        elif line.startswith('A') or line.startswith('C') or line.startswith('G') or line.startswith('T'):
            séquence.replace('/n','')
else:
    
    
        

dictionary = dict(zip(noms, séquence))

#from Bio import SeqIO
#input_file = open("clownfish.fasta")
#my_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

#correction
#open the file
#create a empty dictionary
#iterate over each line
    #line stratswith > --> split

def calcGC(seq):
    gc = (seq.count('G')+seq.count('C'))/len(seq)
    return(gc)
    
def deadFasta(filename):
    fasatDict =[]
    try:
        f=open(filename,'r')
    except:
        print('cannot open the file' %filname)
        
    first=True
    for line in f:
        if line[0] =='>':
            if not first:
                fastaDict[name] = {'seq':dna, 'length':len(dna), 'GC content' : GC}
            tmp =line.split('|')
            desc = tmp[4].split(' ')
            name = desc[1] + ' ' + desc[2] + '_' + tmp[3]
            dna = ''
            first = False
        elif line == /n:
            continue
        else:
            dna = dna + line.rstrip(/n)
    fasta.Dict[name] = {'seq':dna, 'length':len(dna), 'GC content' : GC}
    f.close()
    
#exercies 2

from Bio import SeqIO

input_handle = open("clownfish.gb", "rU")
output_handle = open("clownfish.fasta", "w")

sequences = SeqIO.parse(input_handle, "genbank")
count = SeqIO.write(sequences, output_handle, "fasta")

output_handle.close()
input_handle.close()
print("Converted %i records" % count)


#sans le module:
open('convfasta.fasta','a')
fname = open('clownfish.gb','rU') 
with open('C:\\Users\\karim\\FirstStepProject\\.spyproject\\convfasta.fasta', 'w') as outfile:
    for line in fname:
        with (fname) as infile:
            for line in infile:
                if line.startswith('LOCUS'):
                    outfile.write(line) 
                else:
                    readline()
                if line.startswith('ORIGIN'):
                    while '//' not in line:
                        outfile.write(line)
                else:
                    readline()

