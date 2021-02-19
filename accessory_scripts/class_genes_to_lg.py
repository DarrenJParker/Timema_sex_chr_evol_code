### class_genes_to_lg.py

import sys
import os
import getopt
import decimal
from decimal import *

try:
	opts, args = getopt.getopt(sys.argv[1:], 'l:g:o:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

linkage_file_name = None
gff_file_name = None
out_base = "TEST"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** class_genes_to_lg.py | Written by DJP, 28/11/18 in Python 3.5 in Lausanne  ****\n")
		print("Calcs rel position of genes based on Nosil's linkage map. See 3_RNAseq_mapping.sh")
		
		print("\n***** USAGE *****\n")		
		print("\npython3 class_genes_to_lg.py -l [linkage alignment block file] -g [gff file] -o [output basename]\n\n")
		
		sys.exit(2)
	elif opt in ('-l'):
		linkage_file_name = arg
	elif opt in ('-g'):
		gff_file_name = arg
	elif opt in ('-o'):
		out_base = arg
	else:
		print("i dont know")
		sys.exit(2)
		
		

##############################################

linkage_file = open(linkage_file_name)

linkage_info_dict = {}

seen_scaf = set()

line_N = 0
for line in linkage_file:
	line = line.strip().split("\t")
	line_N = line_N + 1
	if line_N > 1:
		UNIL_scaf = line[0]
		UNIL_coord_start = int(line[3])
		UNIL_coord_end   = int(line[4])
		UNIL_midpoint = int((UNIL_coord_start + UNIL_coord_end) / 2)
		
		
		nosil_coord_start = int(line[9])
		nosil_coord_end   = int(line[10])
		orientation = None
		
		### X chr is not ordered except within scafs - SO make multiple X linkage groups
		lg = line[8]
		if lg == "lgX":
			lg = line[5]
		
		
		
		if UNIL_coord_start > UNIL_coord_end:
			orientation = "reverse"
		else:
			orientation = "forward"
		
		if UNIL_scaf not in seen_scaf:
			seen_scaf.add(UNIL_scaf)
			linkage_info_dict[UNIL_scaf] = [(lg, orientation, nosil_coord_start, nosil_coord_end, UNIL_midpoint)]
		else:
			rec = linkage_info_dict.get(UNIL_scaf)
			rec.append((lg, orientation, nosil_coord_start, nosil_coord_end,UNIL_midpoint ))
			linkage_info_dict[UNIL_scaf] = rec



### select largest alignment block. most have one large block and many little ones.
### should ask kamil about this
### If there was > 1 linkage group I make a note


scaf_to_linkage_info_dict = {}
	

N_single_linkage = 0
N_multi_linkage  = 0

for scaf in linkage_info_dict:
	rec = linkage_info_dict.get(scaf)
	#print(rec)
	lg_set = set()
	lg_set_all = ""

	for i in rec:
		lg_set.add(i[0])
		
	for el in lg_set:
		lg_set_all = lg_set_all + ";" + el
	lg_set_all = lg_set_all.strip(";")

	multi_linkage_groups = None
	if len(lg_set) == 1:
		N_single_linkage = N_single_linkage + 1
		multi_linkage_groups = "NO"
	else:
		N_multi_linkage  = N_multi_linkage + 1
		multi_linkage_groups = "YES"
		
	start_pos = []
	end_pos = []
	curr_block_len = -1
	start_biggest_block_len = None
	end_biggest_block_len   = None
	lg_biggest_block        = None
	orie_biggest_block      = None
	UNIL_mid_biggest_block  = None
	
	for i in rec:
		block_len = i[3] - i[2]
		if curr_block_len < block_len:
			curr_block_len = block_len
			start_biggest_block_len = i[2]
			end_biggest_block_len   = i[3]
			lg_biggest_block   = i[0]
			orie_biggest_block  = i[1]
			UNIL_mid_biggest_block = i[4]
	
	mid_point_biggest_block_len = int((end_biggest_block_len + start_biggest_block_len) / 2)
	
	#print(mid_point_biggest_block_len)	

	scaf_to_linkage_info_dict[scaf] = [mid_point_biggest_block_len,lg_biggest_block,orie_biggest_block,multi_linkage_groups,UNIL_mid_biggest_block,lg_set_all]	


print("\nNumber of scafs assigned to a single linkage group:  " + str(N_single_linkage))	
print("Number of scafs assigned to multiple linkage groups: " + str(N_multi_linkage))	

	

### read gff. work out rel position of each gene.
# NOTE I match the midpoint of the largest alignment block to midpoint of the part of UNIL scaf it aligned to.

out_file = open(out_base + "_gene_rel_pos_lg.csv", "w")
genes_in_lg = 0


out_file.write("gene_name,scaf_name,lg,lg_pos,gene_rel_midpoint,multi_lgs,all_lg\n")

gff_file = open(gff_file_name)
for line in gff_file:
	line = line.strip().split("\t")
	feature = line[2]
	scaf_name = line[0]
	start_gene = int(line[3])
	end_gene = int(line[4])
	gene_mid = int((start_gene + end_gene) / 2)
	
	
	if feature == "gene":
		if scaf_name in seen_scaf:
			
			gene_name = line[8].split("ID=")[1].split(";")[0]
			
			genes_in_lg = genes_in_lg + 1
			linkage_info = scaf_to_linkage_info_dict.get(scaf_name)
			UNIL_mid = linkage_info[4]
			Nosil_mid = linkage_info[0]
			orie = linkage_info[2]
			multi_link = linkage_info[3]
			lg_pos = linkage_info[1]
			all_lg = linkage_info[5]
			lg = lg_pos.split("_")[0]
		
			#print(linkage_info)
			
			
			
			gene_dist_from_UNIL_mid = gene_mid - UNIL_mid 
			# print(line)
			# print(UNIL_mid)
			# print(gene_mid)
			# print(gene_dist_from_UNIL_mid)
			# print(orie)
			# 
			
			gene_mid_rel_to_Nosil = None
			if orie == "forward":
				gene_mid_rel_to_Nosil = Nosil_mid + gene_dist_from_UNIL_mid
			elif orie == "reverse":
				gene_mid_rel_to_Nosil = Nosil_mid - gene_dist_from_UNIL_mid
			else:
				sys.exit(2)
				print("ERROR BAD")
			# print(Nosil_mid)	
			# print(gene_mid_rel_to_Nosil)
			#print(linkage_info)
			
			out_file.write(gene_name + "," + scaf_name + "," + lg + "," + lg_pos + "," + str(gene_mid_rel_to_Nosil) + "," +  multi_link + "," + all_lg + "\n")
			



print("Number of genes in linkage groups: " + str(genes_in_lg))	
	
print("\nFinished, Billy Rock\n\n\n")
	
	
	

