#!/usr/bin/env python
'''---getRNP---'''

import sys
from optparse import OptionParser
import os
import datetime
import subprocess
import itertools
import glob
from collections import deque
from Bio import SeqIO
from Bio.Seq import Seq
from random import randrange

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


def run_rfoot(rfoot_script_path, sam, pred, out):
	print('Running Rfoot')
	if "," in sam: #stranded
		os.system('grep -P "\t\+\t" ' + pred + ' > ' + pred + '_plus')
		os.system('grep -P "\t\-\t" ' + pred + ' > ' + pred + '_minus')
		os.system('perl ' + rfoot_script_path + ' -i ' + sam.split(",")[0] + ' -t ' + pred + '_plus -o tmp/' + out + '_plus_rfoot.txt')
		os.system('perl ' + rfoot_script_path + ' -i ' + sam.split(",")[1] + ' -t ' + pred + '_minus -o tmp/' + out + '_minus_rfoot.txt')
		os.system('cat tmp/' + out + '_plus_rfoot.txt tmp/' + out + '_minus_rfoot.txt > tmp/' + out + '_rfoot.txt')
	else:
		os.system('perl ' + rfoot_script_path + ' -i ' + sam + ' -t ' + pred + ' -o tmp/' + out + '_rfoot.txt')


def separate_reads(sam):
	if "," in sam: #stranded
		os.system('awk \'$1 !~ /@/{ {print > "tmp/n_plus_"length($10)}}\' ' + sam.split(",")[0])
		os.system('awk \'$1 !~ /@/{ {print > "tmp/n_minus_"length($10)}}\' ' + sam.split(",")[1])
	else:
		os.system('awk \'$1 !~ /@/{ {print > "tmp/n_plus_"length($10)}}\' ' + sam)


def read_rfoot(out, thr, thr2):
	out1 = open("tmp/" + out + "_regions.bed","w+")
	regions = {}
	lineas = {}
	n = 0
	for f in ("tmp/" + out + "_plus_rfoot.txt","tmp/" + out + "_minus_rfoot.txt"):
	 	for line in open(f):
	 		if line.startswith("geneID"):
	 			continue
	 		name = line.split("\t")[0]
	 		if not name in regions:
	 			regions[name] = [0,[],line.split("\t")[2],line.split("\t")[3],0]
	 		pme = float(line.split("\t")[11].rstrip("\n"))
	 		n_reads = int(line.split("\t")[8])
	 		if (pme < thr) and (n_reads > thr2):
	 			if line.split("\t")[3] == "+":
	 				out1.write(line.split("\t")[2] + "\t" + str(int(line.split("\t")[4])+1) + "\t" + line.split("\t")[5] + "\t" + name + "_rn" + str(n) + "\t0\t" + line.split("\t")[3] + "\n")
	 				lineas[name + "_rn" + str(n)] = line.split("\t")[2] + "\t" + str(int(line.split("\t")[4])+1) + "\t" + line.split("\t")[5] + "\t" + name + "_rn" + str(n) + "\t0\t" + line.split("\t")[3] + "\n"
	 			else:
	 				out1.write(line.split("\t")[2] + "\t" + line.split("\t")[5] + "\t" + str(int(line.split("\t")[4])+1) + "\t" + name + "_rn" + str(n) + "\t0\t" + line.split("\t")[3] + "\n")
	 				lineas[name + "_rn" + str(n)] = line.split("\t")[2] + "\t" + line.split("\t")[5] + "\t" + str(int(line.split("\t")[4])+1) + "\t" + name + "_rn" + str(n) + "\t0\t" + line.split("\t")[3] + "\n"
	 			n += 1
	out1.close()
	return lineas


def floss_coding(orfs, out, bedtools_path):
	print("FLOSS coding")
	reads = {}
	reads[0] = float(0)
	for file in glob.glob('tmp/*'):
		if file.startswith('tmp/n_') and not file.endswith('.bed') and not file.endswith('_orfs') and not file.endswith('_regions'):
			if not int(file.split("_")[-1]) in reads:
				reads[int(file.split("_")[-1])] = float(0)
			# if "plus" in file:
			# 	os.system('awk \'{print $3"\t"$4"\t"$4+length($10)"\t"$1"\t0\t+"}\' ' + file + '> ' + file.replace("_plus","") + '.bed')
	for file in glob.glob('tmp/*'):
		if file.startswith('tmp/n_') and not file.endswith('.bed') and not file.endswith('_orfs') and not file.endswith('_regions'):
			if not int(file.split("_")[-1]) in reads:
				reads[int(file.split("_")[-1])] = float(0)
			# if "minus" in file:
			# 	os.system('awk \'{print $3"\t"$4"\t"$4+length($10)"\t"$1"\t0\t-"}\' ' + file + '>> ' + file.replace("_minus","") + '.bed')

	# for file in glob.glob('tmp/*'):
	# 	if file.startswith('tmp/n_') and file.endswith('.bed'):
	# 		os.system(bedtools_path + ' coverage -s -b ' + file + ' -a ' + orfs + ' > ' + file + '_orfs')

	for file in glob.glob('tmp/*'):
		if file.endswith('_orfs'):
			n = int(file.replace(".bed_orfs","").split("_")[-1])
			for line in open(file):
				n_reads = float(line.split("\t")[-4])
				reads[n] = reads[n] + n_reads
				reads[0] = reads[0] + n_reads		

	tp = {}
	t = {}
	for n in reads:
		if n == 0:
			continue
		tp[n] = reads[n]/reads[0]
		t[n] = float(0)
		t[0] = float(0)
		print("FLOSS training model: " + str(n) + "\t" + str(tp[n])) 

	return(tp,t)


def floss(tp, t, out, bedtools_path, thr2, thr3, lineas, orfs):
	print("FLOSS main")
	r = {}
	for file in glob.glob('tmp/*'):
		if file.startswith('tmp/n_') and file.endswith('.bed'):
		 	os.system(bedtools_path + ' coverage -s -b ' + file + ' -a tmp/' + out + '_regions.bed > ' + file + '_regions')

	for file in glob.glob('tmp/*'):
		if file.endswith('_regions'):
			n = int(file.replace(".bed_regions","").split("_")[-1])
			for line in open(file):
				n_reads = float(line.split("\t")[-4])
				name = line.split("\t")[3] 
				if not name in r:
					r[name] = {}
					if not n in r[name]:
						r[name][n] = 0
					r[name][n] = r[name][n] + n_reads

	# r = glob.glob('tmp/')
	# for i in r:
	# 	if i.startswith(n_):
	# 		os.remove(i)

	out3 = open("out/" + out + "_regions.bed","w+")
	floss = open("out/" + out + "_regions.floss","w+")
	floss.write("region\tfloss\tn_reads\n")
	done = []
	for name in r:
	 	score = 0
	 	tr = float(sum(r[name]))
	 	if tr < thr2:
	 		continue
	 	for n in tp:
	 		if n != 0:
	 			try:
		 			nr = float(r[name][n])
		 		except:
		 			nr = 0
		 		f = nr / tr
		 		score = score + abs(f - tp[n])
	 	score = score/2
	 	if not name in done:
	 		floss.write(name + "\t" + str(score) + "\t" + str(tr) + "\n")
	 	done.append(name)
	 	if score >= thr3:
	 		out3.write(lineas[name])
	out3.close()
	floss.close()

	os.system(bedtools_path + ' subtract -s -a out/' + out + '_regions.bed -b ' + orfs + ' | sort -k1,1 -k2,2n > out/' + out + '_regions_noorfs.bed')
	os.system(bedtools_path + ' merge -d 0 -s -i out/' + out + '_regions_noorfs.bed -c 4,6 -o distinct,distinct | awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t.\t\"$5}\' > out/' + out + '_regions_merged.bed')


def main():
	usage = "\n%prog  [options]\nNeeded software: Biopython, Rfoot (script in utils), bedtools (2.28 for reproducibility)"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",dest="input",help="Transcript/region file in PRED format (required).")
	parser.add_option("-s","--sam",action="store",dest="sam",help="Ribo-Seq file(s) in SAM format. If stranded data, separate both strands and include both files separated with comma: file+,file-")
	parser.add_option("-c","--cds",action="store",dest="orfs",help="BED file with CDS or translated sequences. Positive model for FLOSS metrics (Required)")
	parser.add_option("-o","--out",action="store",dest="out",help="Output unique name (required).")
	
	rfoot_script = 'utils/Rfoot.pl'
	bedtools_path = 'bedtools'

	(opt,args)=parser.parse_args()

	check_arg(opt.input,"--input")
	check_arg(opt.sam,"--sam")
	check_arg(opt.orfs,"--orfs")
	check_arg(opt.out,"--out")
	check_file(opt.input)
	check_file(opt.sam.split(",")[0])
	check_file(opt.sam.split(",")[1])
	check_file(opt.orfs)

	# run_rfoot(rfoot_script, opt.sam, opt.input, opt.out)
	lineas = read_rfoot(opt.out, 0.6, 10)
	
	# separate_reads(opt.sam)
	(tp,t) = floss_coding(opt.orfs, opt.out, bedtools_path)
	floss(tp, t, opt.out, bedtools_path, 10, 0.35, lineas, opt.orfs)
	
	# r = glob.glob('tmp/')
	# for i in r:
	# 	if i.startswith(out):
	# 		os.remove(i)

if __name__ == '__main__':
	main()

exit(0)