#!/usr/bin/env python
'''---featureCov---'''

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
	usage = "\n%prog  [options]\nNeeded software: Biopython, BLAST (2.3.0 for reproducibility), bedtools (2.28 for reproducibility)"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",dest="regions",help="BED file with region coordinates (output from getRegions.py, required).")
	parser.add_option("-f","--feature",action="store",dest="feature",help="BED/BAM feature to be checked (promoter, CLIP-seq, RNA-seq, Ribo-Seq).")
	parser.add_option("-s","--stranded",action="store",dest="stra",help="Write 'yes' or 'no' for the analysis to require the same strand or not (required).")
	parser.add_option("-o","--out",action="store",dest="out",help="Feature unique name (required).")
	
	bedtools_path = 'bedtools'

	(opt,args)=parser.parse_args()

	check_arg(opt.regions,"--regions")
	check_arg(opt.feature,"--feature")
	check_arg(opt.stra,"--stranded")
	check_arg(opt.out,"--out")
	check_file(opt.regions)
	check_file(opt.feature)

	out2 = open("./tmp/ov.cov","w+")
	if opt.stra == 'yes':
		subprocess.call([bedtools_path, 'coverage', '-s', '-split', '-a', opt.regions, '-b', opt.feature], stdout=out2)
	else:
		subprocess.call([bedtools_path, 'coverage', '-split', '-a', opt.regions, '-b', opt.feature], stdout=out2)
	out2.close()

	#Check coverage and features/kb per region
	regions = {}
	for line in open("./tmp/ov.cov"):
		name = line.split("\t")[3]
		regions.setdefault(name,[0,0,0])
		regions[name][0] = int(line.split("\t")[-4])
		regions[name][1] = int(line.split("\t")[-3])
		regions[name][2] = int(line.split("\t")[-2])

	out = open("out/" + opt.out + ".cov","w+")
	out.write("region\tfeature\tcounts\tlength_cov\tlength\n")
	for name in regions:
		out.write(name + "\t" + opt.out + "\t" + str(regions[name][0]) + "\t" + str(regions[name][1]) + "\t" + str(regions[name][2]) + "\n")
	out.close()


if __name__ == '__main__':
	main()