#!/usr/bin/env python
'''---getRegions---'''

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


class trans_object:
	def __init__(self, chrm, gene, strand, start, end, biotype):
		self.chrm = chrm
		self.gene = gene
		self.strand = strand
		self.start = start
		self.end = end
		self.biotype = biotype


def parse_gtf(gtf, out):
	'''Read a gtf and create a dict with sorted transcript coordinates, chrm, strand, and gene'''
	trans = {}
	bed = open("tmp/" + out + "_total.bed","w+")
	for line in open(gtf):
		if not "\texon\t" in line:
			continue
		t_name = line.split('transcript_id "')[1].split('"')[0]
		g_name = line.split('gene_id "')[1].split('"')[0]

		if 'gene_biotype' in line:
			biot = line.split('gene_biotype "')[1].split('"')[0]
		else:
			biot = "unknown"
		
		trans.setdefault(t_name,trans_object(line.split("\t")[0],g_name,line.split("\t")[6],[],[],biot))
		trans[t_name].start.append(int(line.split("\t")[3]))
		trans[t_name].end.append(int(line.split("\t")[4]))

		bed.write(line.split("\t")[0] + "\t" + line.split("\t")[3] + "\t" + line.split("\t")[4] + "\t" + g_name + "-" + t_name + "\t0\t" + line.split("\t")[6] + "\n")

	[trans[x].start.sort() for x in trans]
	[trans[x].end.sort() for x in trans]	

	bed.close()
	return trans


def getOverlap(intervals):
	sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
	merged = []
	for higher in sorted_by_lower_bound:
		if not merged:
			merged.append(higher)
		else:
			lower = merged[-1]
			if higher[0] <= lower[1] + 30:
				upper_bound = max(lower[1], higher[1])
				merged[-1] = (lower[0], upper_bound)  # replace by merged interval
			else:
				merged.append(higher)
	return merged


def getBiotypes(gtf):
	biotype = {}
	for line in open(gtf):
		if not "\tgene\t" in line:
			continue
		name = line.split('gene_id "')[1].split('"')[0]
		try:
			bt = line.split('gene_biotype "')[1].split('"')[0]
		except:
			bt = "unknown"
		if bt == "protein_coding":
			biotype[name] = "codRNA"
		elif "pseudogene" in bt:
			biotype[name] = "pseudogene"
		elif "IG" in bt:
			biotype[name] = "IG_gene"
		elif "TR" in bt:
			biotype[name] = "TR_gene"
		else:
			biotype[name] = "lncRNA"
	return biotype


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


def run_blast(fasta,db,blast_path,threads,out):
	'''Run BLAST'''
	print("Running BLAST, use version BLAST 2.3.0+ for reproducibility")
	out2 = open("tmp/blast_" + out + ".xml")
	subprocess.call([blast_path, '-query', fasta, '-db', db, '-evalue', '1e-5', '-outfmt', '5', '-num_threads', threads, '-max_target_seqs', '15', '-strand', 'plus', '-word_size', '12'], stdout=out2)
	out2.close()


def parse_blast(code, exc):
	'''Parse BLAST XML output'''
	output_cand = {}
	for line in open(code):
		if '<Iteration_query-def>' in line:
			trans = line.split(">")[1].split("</")[0]
			output_cand[trans] = []
			rev = 0
		elif '<Hit_def>' in line:
			hit_def = line.split(">")[1].split("</")[0]
			aligned_len = 0
		elif '<Hsp_hit-from>' in line:
			st = int(line.split(">")[1].split("</")[0])
		elif '<Hsp_query-frame>' in line:
			if "-" in line.split(">")[1].split("</")[0]:
				rev = 1
		elif '<Hsp_hit-frame>' in line:
			if "-" in line.split(">")[1].split("</")[0]:
				rev = 1			
		elif '<Hsp_hit-to>' in line:
			en = int(line.split(">")[1].split("</")[0])
			if en < st:
				rev = 1
		elif '<Hsp_query-from>' in line:
			start = int(line.split(">")[1].split("</")[0])
		elif '<Hsp_query-to>' in line:
			end = int(line.split(">")[1].split("</")[0])
			if end < start:
				rev = 1
		elif '<Hsp_evalue>' in line:
			e_val = line.split(">")[1].split("</")[0]
		elif '<Hsp_align-len>' in line:
			if int(line.split(">")[1].split("</")[0]) > aligned_len:
				aligned_len = int(line.split(">")[1].split("</")[0])
		
		elif '</Hsp>' in line:
			if trans in exc:
				continue
			if rev == 0:
				#raw_file.write(trans + "\t" + specie + "_homologies_hsa\t" + hit_def.split()[0] + "\t" + gene_t[hit_def.split()[0]] + "\t" + str(aligned_len) + "\t" + str(e_val) + "\n")
				if aligned_len < 30:
					continue
				output_cand[trans].append([aligned_len,e_val,start,end,hit_def.split()[0]])
			rev = 0

	regions = []
	for trans in output_cand:
		match = [[],[],[],1000]
		if len(output_cand[trans]) > 0:
			for elemento2 in output_cand[trans]:
				name = elemento2[4]
				evalue = elemento2[1]
				if float(evalue) < match[3]:
					match[3] = float(evalue)
				s = elemento2[2]
				e = elemento2[3]

				match[1].append((s,e))

		match[2] = getOverlap(match[1])
		for region in match[2]:
			regions.append(trans + ":" + str(region[0]) + "-" + str(region[1]))

	return regions


def exclude(seqs):
	'''Exclude short and masked transcripts'''
	lengths = {}
	exc = []
	for j,elemento in enumerate(seqs):
			if "gene_" in elemento:
					continue

			lengths[elemento] = len(str(seqs[elemento].seq))

			#Filters
			if lengths[elemento] < 200:
					exc.append(elemento)
					continue
			if float(str(seqs[elemento].seq).count("N")) / lengths[elemento] >= 0.75:
					exc.append(elemento)
					continue
			if lengths[elemento] - float(str(seqs[elemento].seq).count("N")) < 100:
					exc.append(elemento)
					continue
	return exc


def write_main_regions(regions, gtf,coords_t, exc, out,mask, ort_list, bedtools_path):
	'''Write a file with all conserved genes according to an external ort_list, f.e. Ensembl Compara. 
	This is only to check overlapping regions with conserved genes; it does not change BLAST output'''
	matches = []
	if ort_list != 'none':
		for line in open(ort_list):
			gene = line.split("\t")[0].rstrip("\n")
			if not gene in matches:
				matches.append(gene)
	
	for region in regions:
		name = region.split(":")[0]
		if not name in matches:
			matches.append(name)

	out_cons = open("tmp/" + out + "_conserved","w+") #File with conserved regions + conserved ort_list of genes
	for line in open(gtf):
		if not "\texon\t" in line:
			continue
		trans = line.split('transcript_id "')[1].split('"')[0]
		if "sense_intronic" in line:
			exc.append(trans)
		if trans in matches:
			out_cons.write(line)
	out_cons.close()

	'''Write main BED file with transcript regions. Afterwards use bedtools to merge and annotate all regions'''
	out2 = open("tmp/" + out + "_regions.bed","w+")
	for region in regions:
		name = region.split(":")[0]
		if name in exc:
			continue
		s = int(region.split(":")[1].split("-")[0])
		e = int(region.split(":")[1].split("-")[1].rstrip("\n"))
		l = sum(coords_t[name].end) - sum(coords_t[name].start) + len(coords_t[name].end)
		if coords_t[name].strand == "-":
			e = l - int(region.split(":")[1].split("-")[0])
			s = l - int(region.split(":")[1].split("-")[1])

		opened = 0
		for n,exon in enumerate(coords_t[name].start):
			e_len = coords_t[name].end[n] - coords_t[name].start[n] + 1	
			if e_len < s:
				s = s - e_len
				e = e - e_len
			else:
				if opened == 0:
					opened = 1
					start = s + coords_t[name].start[n]
				if e_len < e:
					s = s - e_len
					e = e - e_len	
				else:
					end = e + coords_t[name].start[n] - 1
					break

		f = 0
		for n,exon in enumerate(coords_t[name].start):
			if f == 2:
				continue
			elif f == 0:
				if (start >= coords_t[name].start[n]) and (start <= coords_t[name].end[n]):
					if (end >= coords_t[name].start[n]) and (end <= coords_t[name].end[n]):
						out2.write(coords_t[name].chrm + "\t" + str(start) + "\t" + str(end) + "\t" + coords_t[name].gene + "-" + str(region).rstrip("\n") + "\t.\t" + coords_t[name].strand + "\n")
						f = 2
					else:
						out2.write(coords_t[name].chrm + "\t" + str(start) + "\t" + str(coords_t[name].end[n]) + "\t" + coords_t[name].gene + "-" + str(region).rstrip("\n") + "\t.\t" + coords_t[name].strand + "\n")
						f = 1
			elif f == 1:
				if (end >= coords_t[name].start[n]) and (end <= coords_t[name].end[n]):
					out2.write(coords_t[name].chrm + "\t" + str(coords_t[name].start[n]) + "\t" + str(end) + "\t" + coords_t[name].gene + "-" + str(region).rstrip("\n") + "\t.\t" + coords_t[name].strand + "\n")
					f = 2
				else:
					out2.write(coords_t[name].chrm + "\t" + str(coords_t[name].start[n]) + "\t" + str(coords_t[name].end[n]) + "\t" + coords_t[name].gene + "-" + str(region).rstrip("\n") + "\t.\t" + coords_t[name].strand + "\n")
					f = 1
	out2.close()
	
	print("Running bedtools, use version 2.28 for reproducibility")

	os.system('sort -k1,1 -k2,2n tmp/' + out + '_regions.bed > ./tmp/' + out + '_regions.sorted.tmp')
	os.system(bedtools_path + ' merge -d 100 -s -i ./tmp/' + out + '_regions.sorted.tmp -c 4,6 -o distinct,distinct | awk \'{print $1"\t"$2"\t"$3"\t"$4"\t0\t"$5}\' > ./tmp/' + out + '_merged_regions_c.tmp')

	os.system('sort -k1,1 -k2,2n tmp/' + out + '_total.bed > ./tmp/' + out + '_total.sorted.tmp')
	os.system(bedtools_path + ' merge -d 100 -s -i ./tmp/' + out + '_total.sorted.tmp -c 4,6 -o distinct,distinct | awk \'{print $1"\t"$2"\t"$3"\t"$4"\t0\t"$5}\' > tmp/' + out + '_merged_total.bed')
	os.system(bedtools_path + ' subtract -s -a tmp/' + out + '_merged_total.bed -b ./tmp/' + out + '_merged_regions_c.tmp > ./tmp/' + out + '_merged_regions_nc.tmp')

	os.system(bedtools_path + " intersect -S -u -a ./tmp/" + out + "_merged_regions_c.tmp -b tmp/" + out + "_conserved | sed 's/\t0\t/\tov\t/' | sort -k1,1 -k2,2n > ./tmp/" + out + "_merged_regions_ov_r.tmp")
	os.system(bedtools_path + ' merge -d 1 -s -i ./tmp/' + out + '_merged_regions_ov_r.tmp -c 4,6 -o distinct,distinct | awk \'{print $1"\t"$2"\t"$3"\t"$4"\tov\t"$5}\' > ./tmp/' + out + '_merged_regions_ov.tmp') ##

	os.system(bedtools_path + " intersect -S -u -a ./tmp/" + out + "_merged_regions_nc.tmp -b tmp/" + out + "_conserved | sed 's/\t0\t/\tncov\t/' | sort -k1,1 -k2,2n > ./tmp/" + out + "_merged_regions_ovnc_r.tmp")
	os.system(bedtools_path + ' merge -d 1 -s -i ./tmp/' + out + '_merged_regions_ovnc_r.tmp -c 4,6 -o distinct,distinct | awk \'{print $1"\t"$2"\t"$3"\t"$4"\tncov\t"$5}\' > ./tmp/' + out + '_merged_regions_ovnc.tmp')

	os.system(bedtools_path + " subtract -s -a ./tmp/" + out + "_merged_regions_c.tmp -b ./tmp/" + out + "_merged_regions_ov.tmp > ./tmp/" + out + "_merged_regions_c_nov2.tmp")
	
	#Mask non-overlapping regions (f.e. pseudogenes)
	if mask != 'none':
		os.system(bedtools_path + " subtract -a ./tmp/" + out + "_merged_regions_c_nov2.tmp -b " + mask + " > ./tmp/" + out + "_merged_regions_c_nov.tmp")
	else:
		os.system("cp ./tmp/" + out + "_merged_regions_c_nov2.tmp ./tmp/" + out + "_merged_regions_c_nov.tmp")

	os.system(bedtools_path + " subtract -s -a ./tmp/" + out + "_merged_regions_nc.tmp -b ./tmp/" + out + "_merged_regions_ovnc.tmp > ./tmp/" + out + "_merged_regions_nc_nov.tmp")

	os.system("cat ./tmp/" + out + "_merged_regions_c_nov.tmp ./tmp/" + out + "_merged_regions_nc_nov.tmp ./tmp/" + out + "_merged_regions_ov.tmp ./tmp/" + out + "_merged_regions_ovnc.tmp | sort -k1,1 -nk2,2 > ./tmp/"  + out + "_merged_regions_combined.tmp")


def cluster_regions2(out):
	'''CLuster overlapping genes'''
	print('Sortering initial clusters')
	initial_clusters = []
	for line in open("./tmp/" + out + "_merged_regions_combined.tmp"):
		name = line.split("\t")[3]
		c = []
		for n in name.split(","):
			c.append(n.split("-")[0])
		c = list(set(c))
		initial_clusters.append(sorted(set(c)))

	initial_clusters.sort
	lists = list(initial_clusters for initial_clusters,_ in itertools.groupby(initial_clusters))

	resultslist = [] #Create the empty result list.

	print('Merging final clusters')
	if len(lists) >= 1:
		resultlist = [lists[0]]
		if len(lists) > 1:
			for l in lists[1:]:
				listset = set(l)
				merged = False
				for index in range(len(resultlist)):
					rset = set(resultlist[index])
					if len(listset & rset) != 0:
						resultlist[index] = list(listset | rset)
						merged = True
						break
				if not merged: 
					resultlist.append(l)

	clusters = {}
	for result in resultlist:
		for gene in result:
			clusters[gene] = result
	return clusters


def write_final_regions(out,clusters,biotype):
	'''Final output'''
	params = {}
	l = {}
	lines = {}
	for line in open("./tmp/" + out + "_merged_regions_combined.tmp"):
		if (int(line.split("\t")[2]) - int(line.split("\t")[1]) + 1) <= 2:
			continue
		name = line.split("\t")[3]
		clase = "ni"
		if ":" in line:
			if "ov" in line:
				clase = "co"
			else:
				clase = "ci"
		elif "ncov" in line:
			clase = "no"
	
		g = clusters[name.split("-")[0]]

		if not "/".join(g) in params:
			params["/".join(g)] = [0,0,0,0,0,0,0,0]

		if not "/".join(g) in l:
			l["/".join(g)] = 0
			lines["/".join(g)] = []

		l["/".join(g)] = l["/".join(g)] + (int(line.split("\t")[2]) - int(line.split("\t")[1]) + 1)

		bio1v = []
		for e in g:
			if not e.split("-")[0] in biotype:
				bio1v.append("novel")
			else:
				bio1v.append(biotype[e.split("-")[0]])

		if ("pseudogene" in bio1v):
			continue	

		if ("codRNA" in bio1v) or ("IG_gene" in bio1v) or ("TR_gene" in bio1v):
			bio1 = "codRNA"
		elif "lncRNA" in bio1v:
			bio1 = "lncRNA"
		else:
			bio1 = "novel"

		if clase == "ci":
			if params["/".join(g)][0] == 0:
				params["/".join(g)][1] += 1
				params["/".join(g)][0] = 1
				params["/".join(g)][2] = 0
				params["/".join(g)][4] = 0
				params["/".join(g)][6] = 0
			lines["/".join(g)].append(line.split("\t")[0] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\t" + "/".join(g) + "-ci" + str(params["/".join(g)][1]) + "\t" + bio1 + "\t" + line.split("\t")[5])
		elif clase == "co":
			if params["/".join(g)][2] == 0:
				params["/".join(g)][3] += 1
				params["/".join(g)][2] = 1
				params["/".join(g)][0] = 0
				params["/".join(g)][4] = 0
				params["/".join(g)][6] = 0
			lines["/".join(g)].append(line.split("\t")[0] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\t" + "/".join(g) + "-co" + str(params["/".join(g)][3]) + "\t" + bio1 + "\t" + line.split("\t")[5])
		elif clase == "ni":
			if params["/".join(g)][4] == 0:
				params["/".join(g)][5] += 1
				params["/".join(g)][4] = 1
				params["/".join(g)][0] = 0
				params["/".join(g)][2] = 0
				params["/".join(g)][6] = 0
			lines["/".join(g)].append(line.split("\t")[0] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\t" + "/".join(g) + "-ni" + str(params["/".join(g)][5]) + "\t" + bio1 + "\t" + line.split("\t")[5])
		elif clase == "no":
			if params["/".join(g)][6] == 0:
				params["/".join(g)][7] += 1
				params["/".join(g)][6] = 1
				params["/".join(g)][0] = 0
				params["/".join(g)][2] = 0
				params["/".join(g)][4] = 0
			lines["/".join(g)].append(line.split("\t")[0] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\t" + "/".join(g) + "-no" + str(params["/".join(g)][7]) + "\t" + bio1 + "\t" + line.split("\t")[5])

	out2 = open("out/" + out + ".bed","w+")
	for gene in lines:
		if l[gene] >= 30:
			for line in lines[gene]:
				out2.write(line)
	out2.close()


def main():
	usage = "\n%prog  [options]\nNeeded software: Biopython, BLAST (2.3.0 for reproducibility), bedtools (2.28 for reproducibility)"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-f","--fasta",action="store",dest="fasta",help="Transcript file in FASTA format (required).")
	parser.add_option("-g","--gtf",action="store",dest="gtf",help="Transcript coordinates in GTF format (required).")
	parser.add_option("-d","--db",action="store",dest="db",help="FASTA database for species comparison (required).")
	parser.add_option("-o","--out",action="store",dest="out",help="Output unique name (required).")
	parser.add_option("-m","--mask",action="store",dest="mask",default="none",help="GTF file for masking (f.e. pseudogenes).")
	parser.add_option("-O","--orthologs",action="store",dest="ort",default="none",help="Additional list of known gene orthologs (f.e. Ensembl Compara).")
	
	blast_path = 'blastn'
	bedtools_path = 'bedtools'

	(opt,args)=parser.parse_args()

	check_arg(opt.fasta,"--fasta")
	check_arg(opt.gtf,"--gtf")
	check_arg(opt.db,"--db")
	check_arg(opt.out,"--out")
	check_file(opt.fasta)
	check_file(opt.gtf)

	seqs = SeqIO.index(opt.fasta, "fasta") #BLAST 2.3.0
	exc = exclude(seqs)

	#run_blast(opt.fasta,opt.db,blast_path,'1',opt.out)

	print('Parsing')
	regions = parse_blast("tmp/blast_" + opt.out + ".xml", exc)
	coords_t = parse_gtf(opt.gtf,opt.out)

	print('Writing main regions')
	write_main_regions(regions,opt.gtf,coords_t,exc,opt.out,opt.mask,opt.ort,bedtools_path) #BEDtools 2.28
	
	print('Clustering regions from same genes')
	clusters = cluster_regions2(opt.out)

	print('Writing final output')
	biotype = getBiotypes(opt.gtf)
	write_final_regions(opt.out, clusters, biotype)

	r = glob.glob('tmp/*')
	for i in r:
		if i.startswith(out):
			os.remove(i)

if __name__ == '__main__':
	main()

exit(0)
