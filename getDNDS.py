#!/usr/bin/env python
'''---getDNDS---'''

import sys
from optparse import OptionParser
import os
import datetime
import subprocess
import itertools
from collections import deque
from Bio import SeqIO
from Bio.Seq import Seq
from random import randrange
from scipy.stats import chisqprob

__author__ = "Jorge Ruiz-Orera"
__contributor__="Jorge Ruiz-Orera, M.Mar Alba"
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Jorge Ruiz-Orera"
__email__ = "jorruior@gmail.com; malba@imim.es"

def check_arg (arg_str,s):
	'''Check if arg was written'''
	if not arg_str:   # if filename is not given
		print("Error: " + str(s) + " argument not given\n")
		exit()


def check_file (file_str):
	'''Check if input really exists'''
	try:
		open("%s" %file_str)
	except:
		print("Error: " + file_str + " input not found\n")
		exit()


def main():
	usage = "\n%prog  [options]\nNeeded software: Biopython, scipy, PRANK, PAML (v.4)"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-1","--fasta1",action="store",dest="fasta1",help="FASTA file with ORF nucleotides in main species (required).")
	parser.add_option("-2","--fasta2",action="store",dest="fasta2",help="FASTA file with ORF nucleotides in the second species (required).")
	parser.add_option("-o","--out",action="store",dest="out",help="Feature unique name (required).")
	
	prank_folder_path = 'prank'
	paml_folder_path = 'paml'

	(opt,args)=parser.parse_args()

	check_arg(opt.fasta1,"--fasta1")
	check_arg(opt.fasta2,"--fasta2")
	check_arg(opt.out,"--out")
	check_file(opt.fasta1)
	check_file(opt.fasta2)

	seqs1 = SeqIO.index(opt.fasta1, "fasta")
	seqs2 = SeqIO.index(opt.fasta2, "fasta")

	out = open("tmp/" + opt.out + ".aln","w+")
	for orf in seqs1:
		if not orf in seqs2:
			print(orf + " not found in file2")
			continue
		al = [str(seqs1[orf].seq),str(seqs2[orf].seq)]
		out.write(">" + orf + "\n" + al[0] + "\n" + al[1] + "\n")
	out.close()

	out3 = open("out/" + opt.out + "_dnds.tab","w+")
	out4 = open("out/" + opt.out + "_total_homology_final_prot_syn_aligned.aln","w+")
	out3.write("orf\tt\tS\tN\tdnds\tdn\tds\tgaps\tlh\tlhn\tpval\n")
	n = 0
	for line in open("tmp/" + opt.out + ".aln"):
		if ">" in line: 
			s = line.replace(">","").rstrip("\n")
			n = 0
		else:
			n+=1
			if n == 1:
				s1 = line.rstrip("\n")
			elif n == 2:
				s2 = line.rstrip("\n")
				if s2 != "none":
					tmp = open("tmp/tmp.fa","w+")
					tmp.write(">" + s + "_sp1\n" + s1 + "\n>" + s + "_sp2\n" + s2 + "\n")
					tmp.close()

					#Check structure
					p1 = str(Seq(s1).translate(cds=False))
					p2 = str(Seq(s2).translate(cds=False))
					l = len(p2)

					os.system(prank_folder_path + "/src/prank -d=tmp/tmp.fa -o=tmp/tmp_t.fa -translate -F")
					aln1 = ""
					aln2 = ""
					for line in open("tmp/tmp_t.fa.best.nuc.fas"):
						if line.startswith(">"):
							name = line.rstrip("\n")
						elif "_sp1" in name:
							aln1 = aln1 + line.rstrip("\n")	
						elif "_sp2" in name:	
							aln2 = aln2 + line.rstrip("\n")

					while len(aln1) % 3 != 0:
						aln1 = aln1 + "N"
						aln2 = aln2 + "N"

					alnt1 = ""
					alnt2 = ""
					for line in open("tmp/tmp_t.fa.best.pep.fas"):
						if line.startswith(">"):
							name = line.rstrip("\n")
						elif "_sp1" in name:
							alnt1 = alnt1 + line.rstrip("\n")	
						elif "_sp2" in name:	
							alnt2 = alnt2 + line.rstrip("\n")

					c2 = []
					for n,ch in enumerate(alnt1):
						if ch != "-":
							c2.append(n)

					out4.write(">" + s + "\n" + alnt1 + "\n" + alnt2 + "\n")	
					n_gaps = float(alnt1.count("-"))/float(len(alnt1))*100

					#Write FASTA and convert to Phylyp
					fas = open("tmp/fasta.fa","w+")
					fas.write(">" + s + "_sp1\n" + aln1 + "\n")	
					fas.write(">" + s + "_sp2\n" + aln2 + "\n")	
					fas.close()
					os.system("perl " + prank_folder_path + "/Fasta2Phylip.pl tmp/fasta.fa tmp/align.phy")
					os.system("sed -i \'s/\t/  /\' tmp/align.phy")
					
					#Run CODEML
					os.system('echo \"(' + s + '_sp1,' + s + '_sp2)\" > tmp/tree')
					os.system(paml_folder_path + "/bin/codeml utils/codeml.ctl")
					for line in open("tmp/results.codeml"):
						if "lnL" in line:
							lhood = float(line.replace(" ","").split("=")[1].rstrip("\n"))
						elif "dN/dS=" in line:
							line = line.replace("=","= ")
							dnds = line.rstrip("\n").split()
							dnds = [x for x in dnds if not '=' in x]
							break

					os.system(paml_folder_path + "/bin/codeml utils/codeml_fixeddNdS.ctl")
					for line in open("tmp/results_f.codeml"):
						if "lnL" in line:
							lhood_null = float(line.replace(" ","").split("=")[1].rstrip("\n"))
							break
					lhood_ratio = 2*(lhood + abs(lhood_null))
					pval = str(chisqprob(lhood_ratio, 1))			

					if not "-" in dnds:
						out3.write(s + "\t" + "\t".join(map(str,dnds)) + "\t" + str(n_gaps) + "\t" + str(lhood) + "\t" + str(lhood_null) + "\t" + pval + "\n")

	out3.close()
	out4.close()
	os.system("rm 2M*")
	os.system("rm 2N*")
	os.system("rm rst*")
	os.system("rm rub")
	os.system("rm 4fold.nuc")

if __name__ == '__main__':
	main()