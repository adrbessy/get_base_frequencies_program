### ----------------------------
### TFBS Interdistances program
### ----------------------------

'''
This program allows to calculate interdistances between transcription factor binding sites.
You need a matrix with frequency values and fasta sequences (bound sequences (i.e. peaks), and unbound sequences).
This program was written by Adrien Bessy, Arnaud Stigliani and Francois Parcy, and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
#from interdistances_functions import *
import numpy as np
from Bio import SeqIO
import time
import sys 
from operator import truediv
import operator
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
import matplotlib.patches as mpatches
from matplotlib import pylab
import types
import argparse
import logging
from optparse import OptionParser
#from scipy import stats
import collections

parser = argparse.ArgumentParser()                                               

parser.add_argument("--factor", "-fac", type=str, default= "ARF2")
parser.add_argument("--file", "-f", type=str,default= "ARF2_bound_sequences.fas")
parser.add_argument("--matrix", "-mat", type=str,default= "ARF2_score_matrix.txt")
parser.add_argument("--type", "-type", type=str,default= "score")
parser.add_argument("--pseudoCount", "-pc",type=float,default= 0.001)
parser.add_argument("--sequence_number", "-sequence_number", type = int, default= 1000)
args = parser.parse_args()

# python get_base_frequencies.py -fac "ARF2" -pc 0.001 -sequence_number 11654
# python get_base_frequencies.py -fac "ARF5" -pc 0.001 -sequence_number 26659

factorTranscription = args.factor
FastaFile = args.file    
MatrixFile = args.matrix
matrixType = args.type
pseudoCount = args.pseudoCount
sequence_number = args.sequence_number
                    
###################Parameters we can change#################

#factorTranscription = "ARF2_ER7" # ARF2 , ARF5 , LFY_scores_matrix_19nucl , ARF2_ER7 can be choosen

if factorTranscription == "ARF2" :
	FastaFile = "../sequences/ARF2_bound_sequences.fas" 
	MatrixFile = "../matrices/ARF2_matrix1_MEME_1500FirstSeq.txt" #new_ARF2_monomere_matrix_fromER7Matrix.txt, ARF2_OMalley_matrixC.txt
	matrixType = "freq" 
	
if factorTranscription == "ARF5" :
	FastaFile = "../sequences/ARF5_bound_sequences.fas"  
	MatrixFile = "../matrices/ARF5_allSeq_3prime_freq_pasteTo'OMalley.txt" 
	matrixType = "freq" 
#############################################################""

codigoi = { "A" : "T", "C" : "G", "G" : "C", "T" : "A", "N" : "N", "W" : "N"}


def seq_c(site):
        site_i = site[-1::-1]
        site_c = ""
        for x in site_i:
		y = codigoi[x]
                site_c = site_c + y
	return site_c 
  
def get_score_matrix(Mdata,matrixType,pseudoCount):
	## These lines allows to transform the frequency values into scores values
	freq_mat = []
	if matrixType == "count" :
		Mdata = regex.findall(matrix)
		for i in range(0,len(Mdata)):
			if i%4==0:
				#print("Mdata[i] : ",Mdata[i])
				Sum = float(Mdata[i]) + float(Mdata[i+1]) + float(Mdata[i+2]) + float(Mdata[i+3])
				#print("Sum : ",Sum)
				one = float(Mdata[i]) / Sum
				two = float(Mdata[i+1]) / Sum
				three = float(Mdata[i+2]) / Sum
				four = float(Mdata[i+3]) / Sum
				freq_mat.append(one)
				freq_mat.append(two)
				freq_mat.append(three)
				freq_mat.append(four)

		matF = []
		lenMotif=0
		for i in range(0,len(freq_mat)):
			if i%4==0:
				lenMotif=lenMotif+1
				fmax = float(max(freq_mat[i],freq_mat[i+1],freq_mat[i+2],freq_mat[i+3])) + pseudoCount
				for j in range (0,4):
					matF.append(np.log(float(float(freq_mat[i+j]) + pseudoCount) /fmax))
	if matrixType == "freq" :
		Mdata = num.findall(matrix)
		matF = []
		lenMotif=0
		for i in range(0,len(Mdata)):
			if i%4==0:
				lenMotif=lenMotif+1
				fmax = float(max(Mdata[i],Mdata[i+1],Mdata[i+2],Mdata[i+3])) + pseudoCount
				for j in range (0,4):
					matF.append(np.log(float(float(Mdata[i+j]) + pseudoCount) /fmax)) 
	return(matF, lenMotif)
	
def get_DR_basePair_freq(matF,FastaFile,factorTranscription,matRev):
	# This line allows to retrieve all the sequences from the fasta file
	sequences = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))

	print "  There are %s sequence(s) to analyze"%sequence_number
	allScoresPos = []	
	allScoresNeg = []
	sens = ""
	# We will store in these lists all the occurences of each kind of interdistances between motifs found in all sequences.
	index = 0
	#bestScoreBySeq1 = []
	#bestScoreBySeq2 = []
	#bestScoreBySeq3 = []
	#bestScoreBySeq4 = []
	#bestScoreBySeq5 = []
	#bestScoreBySeq6 = []
	#bestScoreBySeq7 = []
	#bestScoreBySeq8 = []
	#bestScoreBySeq9 = []
	#bestScoreBySeq10 = []
	#bestScoreBySeq11 = []
	#bestScoreBySeq12 = []
	bestScoreBySeq13 = []
	#bestScoreBySeq14 = []
	#bestScoreBySeq15 = []
	#bestScoreBySeq16 = []
	nb = 0
	# We look at all the fasta sequences:
	for s in sequences:
			# We will store in this list all the best scores (see the threshold after) found for subsequences of one sequence
			#if type(threshold) is list:
			bestScore = 0
			positionOfTheBestScore = 0
			# This line allows to retrieve the DNA sequence
			seq = sequences[s].seq
			id=sequences[s].id
			#score_seq1 = []
			#score_seq2 = []
			#score_seq3 = []
			#score_seq4 = []
			#score_seq5 = []
			#score_seq6 = []
			#score_seq7 = []
			#score_seq8 = []
			#score_seq9 = []
			#score_seq10 = []
			#score_seq11 = []
			#score_seq12 = []
			score_seq13 = []
			#score_seq14 = []
			#score_seq15 = []
			#score_seq16 = []
			# We look at each sub-sequences of the whole sequence. Each sub-sequence has the same length that the matrix length.
			for c in range(len(seq) - (lenMotif +10)):
				strandPos = seq[c:c+lenMotif].upper()
				#print("strandPos : ",strandPos)
				strandPos_ext = seq[c-10:c+lenMotif].upper()
				#print("strandPos : ",strandPos)
				#strandNeg = seq_c(strandPos)
				strandNeg_ext = seq[c:c+lenMotif+11].upper()
				#print("strandNeg : ",strandNeg)
				test = 0
				for nu in strandPos :
					if nu not in ["A","C","G","T"]:
						test = 1
				for nu in strandPos_ext :
					if nu not in ["A","C","G","T"]:
						test = 1
				for nu in strandNeg_ext :
					if nu not in ["A","C","G","T"]:
						test = 1
				for nu in seq_c(seq[c-1].upper()) :
					if nu not in ["A","C","G","T"]:
						test = 1
				if test == 1:
					score = "NA"
				else :
					n = 0
					#These lines allows to calculate a score for one sub-sequence
					scoreStrandPos = 0
					scoreStrandNeg = 0
					while n<lenMotif:
						if strandPos[n] == 'A':
							scoreStrandPos = scoreStrandPos + matF[n*4]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4]
						elif strandPos[n] == 'C':
							scoreStrandPos = scoreStrandPos + matF[n*4+1]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4+1]
						elif strandPos[n] == 'G':
							scoreStrandPos = scoreStrandPos + matF[n*4+2]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4+2]
						elif strandPos[n] == 'T':
							scoreStrandPos = scoreStrandPos + matF[n*4+3]
							scoreStrandNeg = scoreStrandNeg + matRev[n*4+3]
						n += 1
					#score_seq1.append([scoreStrandPos,str(seq[c-10].upper())])
					#score_seq1.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif+10].upper()))])
					#score_seq2.append([scoreStrandPos,str(seq[c-9].upper())])
					#score_seq2.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif+9].upper()))])
					#score_seq3.append([scoreStrandPos,str(seq[c-8].upper())])
					#score_seq3.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif+8].upper()))])
					#score_seq4.append([scoreStrandPos,str(seq[c-7].upper())])
					#score_seq4.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif+7].upper()))])
					#score_seq5.append([scoreStrandPos,str(seq[c-6].upper())])
					#score_seq5.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif+6].upper()))])
					#score_seq6.append([scoreStrandPos,str(seq[c-5].upper())])
					#score_seq6.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif+5].upper()))])
					#score_seq7.append([scoreStrandPos,str(seq[c-4].upper())])
					#score_seq7.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif+4].upper()))])
					#score_seq8.append([scoreStrandPos,str(seq[c-3].upper())])
					#score_seq8.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif+3].upper()))])
					#score_seq9.append([scoreStrandPos,str(seq[c-2].upper())])
					#score_seq9.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif+2].upper()))])
					#score_seq10.append([scoreStrandPos,str(seq[c-1].upper())])
					#score_seq10.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif+1].upper()))])
					#score_seq11.append([scoreStrandPos,str(seq[c].upper())])
					#score_seq11.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif].upper()))])
					#score_seq12.append([scoreStrandPos,str(seq[c+1].upper())])
					#score_seq12.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif-1].upper()))])
					#score_seq13.append([scoreStrandPos,str(seq[c+2].upper())])
					#score_seq13.append([scoreStrandNeg,str(seq_c(seq[c+lenMotif-3].upper()))])
					#score_seq13.append([scoreStrandPos,str(seq[c+5:c+7].upper())])
					#score_seq13.append([scoreStrandNeg,str(seq_c(seq[c+4:c+6].upper()))])
					score_seq13.append([scoreStrandPos,str(seq[c+7:c+9].upper())])
					score_seq13.append([scoreStrandNeg,str(seq_c(seq[c+2:c+4].upper()))])
					#score_seq14.append([scoreStrandPos,str(seq[c+lenMotif+1].upper())])
					#score_seq14.append([scoreStrandNeg,str(seq_c(seq[c-1].upper()))])
					#score_seq15.append([scoreStrandPos,str(seq[c+lenMotif+2].upper())])
					#score_seq15.append([scoreStrandNeg,str(seq_c(seq[c-2].upper()))])
					#score_seq16.append([scoreStrandPos,str(seq[c+lenMotif+3].upper())])
					#score_seq16.append([scoreStrandNeg,str(seq_c(seq[c-3].upper()))])
			#bestScoreBySeq1.append(max(score_seq1, key=lambda x: x[0]))
			#bestScoreBySeq2.append(max(score_seq2, key=lambda x: x[0]))
			#bestScoreBySeq3.append(max(score_seq3, key=lambda x: x[0]))
			#bestScoreBySeq4.append(max(score_seq4, key=lambda x: x[0]))
			#bestScoreBySeq5.append(max(score_seq5, key=lambda x: x[0]))
			#bestScoreBySeq6.append(max(score_seq6, key=lambda x: x[0]))
			#bestScoreBySeq7.append(max(score_seq7, key=lambda x: x[0]))
			#bestScoreBySeq8.append(max(score_seq8, key=lambda x: x[0]))
			#bestScoreBySeq9.append(max(score_seq9, key=lambda x: x[0]))
			#bestScoreBySeq10.append(max(score_seq10, key=lambda x: x[0]))
			#bestScoreBySeq11.append(max(score_seq11, key=lambda x: x[0]))
			#bestScoreBySeq12.append(max(score_seq12, key=lambda x: x[0]))
			bestScoreBySeq13.append(max(score_seq13, key=lambda x: x[0]))
			#bestScoreBySeq14.append(max(score_seq14, key=lambda x: x[0]))
			#bestScoreBySeq15.append(max(score_seq15, key=lambda x: x[0]))
			#bestScoreBySeq16.append(max(score_seq16, key=lambda x: x[0]))

			index = index + 1
			if sequence_number :
				nb = nb + 1
			if nb == sequence_number : 
				break
	#bs1 = [item[1] for item in bestScoreBySeq1]
	#bs2 = [item[1] for item in bestScoreBySeq2]
	#bs3 = [item[1] for item in bestScoreBySeq3]
	#bs4 = [item[1] for item in bestScoreBySeq4]
	#bs5 = [item[1] for item in bestScoreBySeq5]
	#bs6 = [item[1] for item in bestScoreBySeq6]
	#bs7 = [item[1] for item in bestScoreBySeq7]
	#bs8 = [item[1] for item in bestScoreBySeq8]
	#bs9 = [item[1] for item in bestScoreBySeq9]
	#bs10 = [item[1] for item in bestScoreBySeq10]
	#bs11 = [item[1] for item in bestScoreBySeq11]
	#bs12 = [item[1] for item in bestScoreBySeq12]
	bs13 = [item[1] for item in bestScoreBySeq13]
	#bs14 = [item[1] for item in bestScoreBySeq14]
	#bs15 = [item[1] for item in bestScoreBySeq15]
	#bs16 = [item[1] for item in bestScoreBySeq16]
	return(bs13)
	#return(bs1,bs2,bs3,bs4,bs5,bs6,bs7,bs8,bs9,bs10,bs11,bs12,bs13)

########################################### About the main matrix #######################

''' The sens of the matrix is important: The positions are on the vertical sens and the bases are on the horizontal sens as described in the example.
separation between numbers can be spaces, tabulation, comas...

                                                                         Example :   A C G T
                                                                  position 1:           0.16456   0.21614       0.1565,0.1645
                                                                  position 2:           0.645; 0.654    0.155 |||||| 0.4444
                                                                                        ...
                                                                        '''
####################################################################################

# These 3 lines allows to retrieve the matrix from the file
F = open(MatrixFile,"r")
matrix = F.read().replace("\r","\n") + "\n"
F.close()

# These 3 lines allows to retrieve all the individual frequency values from the matrix and put them in order into a list
import re
num = re.compile(r"([+-]?\d+[.,]\d+)")
Mdata = num.findall(matrix)

matScore, lenMotif = get_score_matrix(Mdata,matrixType,pseudoCount)

# The following line allows to produce the reversed matrix
'''if we take the example given before : A T G C
			Position 1:      0.4444  0.155  0.654   0.645
			Position 2:      0.1645  0.1565 0.21614 0.16456
Now, we can notice that scores change between the positions 1 and 2, and between A and T, and between G and C.
So we can calculate with this reverse matrix, the score of the complementary strand.
'''
matRev = list(reversed(matScore))

########## get INTERDISTANCE VALUES for POSITIVE sets:

#bs1,bs2,bs3,bs4,bs5,bs6,bs7,bs8,bs9,bs10,bs11,bs12,bs13 = get_DR_basePair_freq(matScore,FastaFile,threshold,factorTranscription,Interdistance_maxValue,matRev)
bs13 = get_DR_basePair_freq(matScore,FastaFile,factorTranscription,matRev)

#c1 =  {n: float(bs1.count(n))/float(len(bs1)) for n in bs1}
##c1 = sorted(c1.items(), key=operator.itemgetter(1),reverse=True)
#c1 = collections.OrderedDict(sorted(c1.items()))
#print("c1 : ",c1,"\n")

#c2 =  {n: float(bs2.count(n))/float(len(bs2)) for n in bs2}
##c2 = sorted(c2.items(), key=operator.itemgetter(1),reverse=True)
#c2 = collections.OrderedDict(sorted(c2.items()))
#print("c2 : ",c2,"\n")

#c3 =  {n: float(bs3.count(n))/float(len(bs3)) for n in bs3}
##c3 = sorted(c3.items(), key=operator.itemgetter(1),reverse=True)
#c3 = collections.OrderedDict(sorted(c3.items()))
#print("c3 : ",c3,"\n")

#c4 =  {n: float(bs4.count(n))/float(len(bs4)) for n in bs4}
##c4 = sorted(c4.items(), key=operator.itemgetter(1),reverse=True)
#c4 = collections.OrderedDict(sorted(c4.items()))
#print("c4 : ",c4,"\n")

#c5 =  {n: float(bs5.count(n))/float(len(bs5)) for n in bs5}
##c5 = sorted(c5.items(), key=operator.itemgetter(1),reverse=True)
#c5 = collections.OrderedDict(sorted(c5.items()))
#print("c5 : ",c5,"\n")

#c6 =  {n: float(bs6.count(n))/float(len(bs6)) for n in bs6}
##c6 = sorted(c6.items(), key=operator.itemgetter(1),reverse=True)
#c6 = collections.OrderedDict(sorted(c6.items()))
#print("c6 : ",c6,"\n")

#c7 =  {n: float(bs7.count(n))/float(len(bs7)) for n in bs7}
##c7 = sorted(c7.items(), key=operator.itemgetter(1),reverse=True)
#c7 = collections.OrderedDict(sorted(c7.items()))
#print("c7 : ",c7,"\n")

#c8 =  {n: float(bs8.count(n))/float(len(bs8)) for n in bs8}
##c8 = sorted(c8.items(), key=operator.itemgetter(1),reverse=True)
#c8 = collections.OrderedDict(sorted(c8.items()))
#print("c8 : ",c8,"\n")

#c9 =  {n: float(bs9.count(n))/float(len(bs9)) for n in bs9}
##c9 = sorted(c9.items(), key=operator.itemgetter(1),reverse=True)
#c9 = collections.OrderedDict(sorted(c9.items()))
#print("c9 : ",c9,"\n")

#c10 =  {n: float(bs10.count(n))/float(len(bs10)) for n in bs10}
#c10 = collections.OrderedDict(sorted(c10.items()))
#print("c10 : ",c10,"\n")

#c11 =  {n: float(bs11.count(n))/float(len(bs11)) for n in bs11}
##c11 = sorted(c11.items(), key=operator.itemgetter(1),reverse=True)
#c11 = collections.OrderedDict(sorted(c11.items()))
#print("c11 : ",c11,"\n")

#c12 =  {n: float(bs12.count(n))/float(len(bs12)) for n in bs12}
##c12 = sorted(c12.items(), key=operator.itemgetter(1),reverse=True)
#c12 = collections.OrderedDict(sorted(c12.items()))
#print("c12 : ",c12,"\n")

c13 =  {n: float(bs13.count(n))/float(len(bs13)) for n in bs13}
c13 = sorted(c13.items(), key=operator.itemgetter(1),reverse=True)
#c13 = collections.OrderedDict(sorted(c13.items()))
print("c13 : ",c13,"\n")

#c14 =  {n: float(bs14.count(n))/float(len(bs14)) for n in bs14}
#c14 = collections.OrderedDict(sorted(c14.items()))
#print("c14 : ",c14,"\n")

#c15 =  {n: float(bs15.count(n))/float(len(bs15)) for n in bs15}
#c15 = collections.OrderedDict(sorted(c15.items()))
#print("c15 : ",c15,"\n")

#c16 =  {n: float(bs16.count(n))/float(len(bs16)) for n in bs16}
#c16 = collections.OrderedDict(sorted(c16.items()))
#print("c15 : ",c15,"\n")

#content = "A\tC\tG\tT\n\n"+str(c1)+"\n"+str(c2)+"\n"+str(c3)+"\n"+str(c4)+"\n"+str(c5)+"\n"+str(c6)+"\n"+str(c7)+"\n"+str(c8)+"\n"+str(c9)+"\n"+str(c10)+"\n"+str(c11)+"\n"+str(c12)+"\n"+str(c13)
content = "A\tC\tG\tT\n\n"+str(c13)

text_file = open(factorTranscription+"_twoLastBasesFreq.txt", "w")
text_file.write(str(content).replace("OrderedDict([('A', ","").replace(")])","").replace("), ('C', ","\t").replace("), ('G', ","\t").replace("), ('T', ","\t"))
text_file.close()			
		
