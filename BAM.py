#!/usr/bin/env python
# coding: utf-8

"""ORF FINDER. This programm takes in fasta-file, executes ORF-search and writes all ORFs, that have been found into CSV-file"""

"""ORF (open reading frame) ist ein DNA Abschnitt zwischen Start- und Stoppcodon. ORF-Finder gibt alle möglichen ORFs heraus. 
Bei der prokariotischen DNA-Sequenz sind insgesamt 6 ORFs möglich. Jedoch wird nur eine von der 6 möglichen Reading Frames ein ORF.
Normalerweise ist die längste reading Frame ein ORF."""

import re
import Bio
from Bio import SeqIO
import csv

"""DATA PREPARATION"""

# Als erstes muss FASTA File geparst werden: "/home/maria/.local/share/applications/gencode.v41.pc_transcripts.fa"
name = input("path to FASTA-file: ")
sequences = SeqIO.parse(name,"fasta") # iterierbares Objekt
print("your FASTA-file: ",sequences)

# Looping through all records in FASTA-file
for record in sequences:
    transcript_1 = record
    break
print(transcript_1) # die Attributen von diesem bestimmten Objekt werden wiedergegeben
# die erhaltene Sequenz ist in Ein-Buchstabencode eingespeichert. Nun muss es in Zero-basierende Form umgeschrieben werden.

# Diese seq enthält UTRs. Man braucht nur CDS.
CDS = transcript_1.seq[60:1041]
print("CDS: ", CDS, type(CDS)) # es beginnt mit ATG end endet mit TAG.

# type of the CDS object ist Class. We need is as a string. Before that we create a list with all records separted by ">"-sign.
L = []
for record in sequences:
    L.append(record)
print("length: ",len(L), "type: ",type(L), "1st SeqRecord (list element): ",L[0], "type of SeqRecord:", type(L[0]))
# nun werden die SeqRecord Objects durch Indexing vom List wiedergegeeben.

# Typ vom List Object an der Position 0 ist SeqRecord Object. Nun kann die Sequenz durch Methods abgerufen werden.
print("string:",L[0].seq, type(L[0].seq))

seq_as_string= str(L[0].seq)
my_seq = seq_as_string

"""FIND ORFs"""

#defining a function which
def find_ORFs(my_seq):
    ORFs = []
    if "ATG" in my_seq:
        for startMatch in re.finditer("ATG", my_seq):
            remaining = my_seq[startMatch.start():]
            #print(remaining)
            for stopMatch in re.finditer("TAA|TGA|TAG", remaining):
                substring = remaining[:stopMatch.end()]
                #print(substring)
                if len(substring) % 3 == 0:
                    #print(substring)
                    ORFs.append(substring)
                    break
    #print(ORFs)
    ORFs.sort(key=len, reverse=True)
    #print(ORFs)
    #return ORFs[0] returns the longes orf
    return ORFs

result = find_ORFs(my_seq)
print("ORFs that has been found: ", result)
data = result

# nun kann das entsprechende File im PC-Ordner wiedergefunden werden
with open("results.csv", "w") as file:
    writer = csv.writer(file)
    writer.writerow(data)




