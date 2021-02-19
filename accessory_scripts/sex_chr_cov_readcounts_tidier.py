### sex_chr_cov_readcounts_tidier.py

import sys
import os
import getopt
import decimal
from decimal import *

try:
	opts, args = getopt.getopt(sys.argv[1:], 'o:r:l:y:c:p:d:g:L:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


out_prefix = "Testttttt"
readcount_filename = None
gene_len_filename = None
M_F_cov_filename = None
gff_filename  = None
orth_matrix_filename = None
linkage_filename = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** sex_chr_cov_readcounts_tidier.py | Written by DJP, 30/07/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Adds chromosome info etc to read count files\nSee 3_RNAseq_mapping.sh for actual usage")
		
		print("\n**** USAGE **** \n")
		print("python3 sex_chr_cov_readcounts_tidier.py -r [read_count file] -l [gene length file] -c [M/F coverage file] -y [orth_matrix_filename] -L [linkage_filename] -g [gff file] -o [out prefix] [options] \n")
		
		print("\n**** OPTIONS ****\n")
		print("-o\toutput prefix")
		print("-r\tread_count file - output from HTSeq_to_edgeR.py")
		print("-l\tgene lengths (Total exon length by gene) - output from gff_feature_lengths.py")
		print("-c\tMale/female coverage file - output from Male_and_female_coverage.R")
		print("-y\tOtholog matrix file")
		print("-L\tLinkage info file - output from class_genes_to_lg.py ")	
		print("-g\tgff file - Use the one outputted by Maker_gff_to_HTseq_gff.py\n\n")

		sys.exit(2)

	elif opt in ('-o'):
		out_prefix = arg
	elif opt in ('-r'):
		readcount_filename = arg
	elif opt in ('-l'):
		gene_len_filename  = arg
	elif opt in ('-y'):
		orth_matrix_filename  = arg
	elif opt in ('-c'):
		M_F_cov_filename = arg
	elif opt in ('-g'):
		gff_filename = arg
	elif opt in ('-L'):
		linkage_filename = arg
	else:
		print("i dont know")
		sys.exit(2)


if readcount_filename == None:
	print("\nError: missing readcount_filename")
	sys.exit(2)
	
if gene_len_filename == None:
	print("\nError: missing gene_len_filename")
	sys.exit(2)
	
if M_F_cov_filename == None:
	print("\nError: missing M_F_cov_filename")
	sys.exit(2)
	
if gff_filename == None:
	print("\nError: missing gff_filename")
	sys.exit(2)
	
if orth_matrix_filename== None:
	print("\nError: missing orth_matrix_filename")
	sys.exit(2)

if linkage_filename == None:
	print("\nError: missing linkage_filename")
	sys.exit(2)


out_file_stats = open(out_prefix + "_chr_info_and_counts_RUNINFO.txt", "w")
out_file_stats.write("Files:\n" + str(readcount_filename) + "\n") 
out_file_stats.write(str(gene_len_filename) + "\n")
out_file_stats.write(str(M_F_cov_filename) + "\n")
out_file_stats.write(str(gff_filename) + "\n")
out_file_stats.write(str(linkage_filename) + "\n")
out_file_stats.write(str(orth_matrix_filename) + "\n\n")

################################################################################################################################################
######## read M/F coverage file

M_F_cov_file = open(M_F_cov_filename)

M_F_dict = {}

X_scafs_N_hard = 0
X_len_hard     = 0
A_scafs_N_hard = 0
A_len_hard     = 0

X_scafs_N_soft = 0
X_len_soft     = 0
A_scafs_N_soft = 0
A_len_soft     = 0

line_N = 0
for line in M_F_cov_file:
	line_N = line_N + 1
	line = line.rstrip("\n").replace('"', '').split(",")
	if line_N  == 1:
		contig_name_index = line.index("contig")
		contig_length_index = line.index("length")
		M_F_mode_index = line.index("M_F_mode")
		class_soft_index = line.index("class_soft")
		class_hard_index = line.index("class_hard")
		
		#print(line)
	
	if line_N > 1:
		contig_name   = line[contig_name_index]
		contig_length = line[contig_length_index]
		M_F_mode      = line[M_F_mode_index]
		class_soft    = line[class_soft_index]
		class_hard    = line[class_hard_index]

		if class_soft == "X":
			X_scafs_N_soft = X_scafs_N_soft + 1
			X_len_soft     = X_len_soft + int(contig_length)
		
		elif class_soft == "A":
			A_scafs_N_soft = A_scafs_N_soft + 1
			A_len_soft     = A_len_soft + int(contig_length)
			
		else:
			print("Error")
			sys.exit(2)
		
		if class_hard == "X":
			X_scafs_N_hard = X_scafs_N_hard + 1
			X_len_hard     = X_len_hard + int(contig_length)
		
		elif class_hard == "A":
			A_scafs_N_hard = A_scafs_N_hard + 1
			A_len_hard     = A_len_hard + int(contig_length)
			
		else:
			print("Error")
			sys.exit(2)
					
		M_F_dict[contig_name] = [contig_length, M_F_mode, class_soft, class_hard]

M_F_cov_file.close()



print("\n\nN Autosomal scafs (soft class): " + str(A_scafs_N_soft))
print("Number of bases Autosomal (soft class):  " + str(A_len_soft))
print("N X-linked scafs (soft class): " + str(X_scafs_N_soft))
print("Number of bases X-linked (soft class):  " + str(X_len_soft))
print("Percentage of genome X-linked  (soft class): " + str((X_len_soft / (X_len_soft + A_len_soft)) * 100) + "\n")

print("N Autosomal scafs (hard class): " + str(A_scafs_N_hard))
print("Number of bases Autosomal (hard class):  " + str(A_len_hard))
print("N X-linked scafs (hard class): " + str(X_scafs_N_hard))
print("Number of bases X-linked (hard class):  " + str(X_len_hard))
print("Percentage of genome X-linked  (hard class): " + str((X_len_hard / (X_len_hard + A_len_hard)) * 100) + "\n")

out_file_stats.write("N Autosomal scafs (soft class): " + str(A_scafs_N_soft) + "\n")
out_file_stats.write("Number of bases Autosomal (soft class):  " + str(A_len_soft) + "\n")
out_file_stats.write("N X-linked scafs (soft class): " + str(X_scafs_N_soft) + "\n")
out_file_stats.write("Number of bases X-linked (soft class):  " + str(X_len_soft) + "\n")
out_file_stats.write("Percentage of genome X-linked  (soft class): " + str((X_len_soft / (X_len_soft + A_len_soft)) * 100) + "\n")

out_file_stats.write("N Autosomal scafs (hard class): " + str(A_scafs_N_hard) + "\n")
out_file_stats.write("Number of bases Autosomal (hard class):  " + str(A_len_hard) + "\n")
out_file_stats.write("N X-linked scafs (hard class): " + str(X_scafs_N_hard) + "\n")
out_file_stats.write("Number of bases X-linked (hard class):  " + str(X_len_hard) + "\n")
out_file_stats.write("Percentage of genome X-linked  (hard class): " + str((X_len_hard / (X_len_hard + A_len_hard)) * 100) + "\n")
 
################################################################################################################################################
######## gene name to scaff dict

# using gff

gene_to_scaf_dict = {}

gff_file = open(gff_filename)
for line in gff_file:
	line = line.rstrip("\n")
	feature = line.split("\t")[2]
	scaf_n = line.split("\t")[0]
	if feature == "exon":
		descrip_line = line.split("\t")[8]
		gene_id = descrip_line.split("gene_id=")[1].split(";")[0]
		gene_to_scaf_dict[gene_id] = scaf_n
		
gff_file.close()


###############################################################################################################################################
###### read gene lengths into dict

gene_len_file = open(gene_len_filename)

gene_len_dict = {}

line_N = 0
for line in gene_len_file:
	line_N = line_N + 1
	if line_N > 1:
		line = line.rstrip("\n").split(",")
		#print(line)
		gene_len_dict[line[0]] = line[1]
gene_len_file.close()


###############################################################################################################################################
###### linkage info


linkage_file = open(linkage_filename)

LG_dict = {}

line_N = 0
for line in linkage_file:
	line_N = line_N + 1
	if line_N > 1:
		line = line.rstrip("\n").split(",")
		all_lg = line[6].split(";")
		lgs = set()
		for l in all_lg:
			l = l.split("_")[0]
			lgs.add(l)
		
		### if different scafs from lgX - class as multi
		multi_lg = "NA"
		if len(lgs) == 1:
			multi_lg = "NO"
		else:
			multi_lg = "YES"

		LG_dict[line[0]] = [line[3], line[4], multi_lg, line[6]]

linkage_file.close()


######################################################################################################
### orth info

gene_to_orth_dict = {}
orth_matrix_file = open(orth_matrix_filename)
for line in orth_matrix_file:
	line = line.rstrip("\n").split("\t")
	HOG_name = line[0]
	
	for i in range(1, len(line)):
		gene_name = line[i]
		gene_to_orth_dict[gene_name] = HOG_name

	
orth_matrix_file.close()

################################################################################################################################################
######## read readcount file and export

readcount_file = open(readcount_filename)
out_file = open(out_prefix + "_chr_info_and_counts.csv", "w")

line_N = 0
gene_no_scaf_info = 0
genes_X_hard = 0
genes_A_hard = 0
genes_X_soft = 0
genes_A_soft = 0

N_orths = 0

genes_X_soft_orth = 0
genes_A_soft_orth = 0
genes_X_hard_orth = 0
genes_A_hard_orth = 0


for line in readcount_file:
	line_N = line_N + 1
	line = line.rstrip("\n").split(",")
	
	counts = ""
	for i in range(1,len(line)):
		counts = counts + "," + line[i]
	
	gene_id = line[0]
	
	if line_N == 1:
		out_file.write("gene_id,tot_exon_len,HOG,scaf_name,scaf_len,M_F_mode,chr_soft,chr_hard" + counts + ",lg_pos,gene_rel_midpoint,multi_lg,all_lg\n")
		
	else:
		
		scaf_name = gene_to_scaf_dict.get(gene_id)
		if scaf_name == None:
			print("\ncant find scaf name for some genes. Exiting\n")
			sys.exit(2)

		scaf_len = ""
		M_F_mode = ""
		chr_hard = ""
		chr_soft = ""
		
		scaf_rec = M_F_dict.get(scaf_name)
		if scaf_rec == None:
			scaf_len = "NA"
			M_F_mode = "NA"
			chr_hard = "NA"
			chr_soft = "NA"
			gene_no_scaf_info = gene_no_scaf_info + 1
		else:
			scaf_len = scaf_rec[0]
			M_F_mode = scaf_rec[1]
			chr_soft = scaf_rec[2]
			chr_hard = scaf_rec[3]

			if chr_soft == "X":
				genes_X_soft = genes_X_soft + 1
			if chr_soft == "A":
				genes_A_soft = genes_A_soft + 1

			if chr_hard == "X":
				genes_X_hard = genes_X_hard + 1
			if chr_hard == "A":
				genes_A_hard = genes_A_hard + 1
		
		## gene len			
		gene_len = gene_len_dict.get(gene_id)
		if gene_len == None:
			print("\ncant find gene length for some genes. Exiting\n")
			sys.exit(2)
		
		## lg info
		lg_rec = LG_dict.get(gene_id)
		if lg_rec == None:
			lg_rec = ["NA", "NA", "NA", "NA"]
		
		lg_out = ""
		for l in lg_rec:
			lg_out = lg_out + "," + l

		### orth info
		
		HOG_name = gene_to_orth_dict.get(gene_id)
		if HOG_name == None:
			HOG_name = "NA"
		else:
			N_orths = N_orths + 1
			if chr_soft == "X":
				genes_X_soft_orth = genes_X_soft_orth + 1
			if chr_soft == "A":
				genes_A_soft_orth = genes_A_soft_orth + 1

			if chr_hard == "X":
				genes_X_hard_orth = genes_X_hard_orth + 1
			if chr_hard == "A":
				genes_A_hard_orth = genes_A_hard_orth + 1			
		
		out_file.write(gene_id + "," + str(gene_len) + "," +  HOG_name + "," + scaf_name + "," + str(scaf_len) + "," + M_F_mode +  "," + chr_soft + "," + chr_hard + counts + lg_out + "\n")


readcount_file.close()


print("\nNumber of genes Autosomal (soft): " + str(genes_A_soft))
print("Number of genes X-linked (soft): " + str(genes_X_soft))
print("Number of genes Autosomal (hard): " + str(genes_A_hard))
print("Number of genes X-linked (hard): " + str(genes_X_hard))
print("Number of genes with no chromosome info: " + str(gene_no_scaf_info))

print("\nNumber of orths: " + str(N_orths ))
print("Number of orthologs Autosomal (soft): " + str(genes_A_soft_orth))
print("Number of orthologs  X-linked (soft): " + str(genes_X_soft_orth))
print("Number of orthologs  Autosomal (hard): " + str(genes_A_hard_orth))
print("Number of orthologs  X-linked (hard): " + str(genes_X_hard_orth))

out_file_stats.write("\nNumber of genes Autosomal (soft): " + str(genes_A_soft)+ "\n")
out_file_stats.write("Number of genes X-linked (soft): " + str(genes_X_soft)+ "\n")
out_file_stats.write("Number of genes Autosomal (hard): " + str(genes_A_hard)+ "\n")
out_file_stats.write("Number of genes X-linked (hard): " + str(genes_X_hard)+ "\n")
out_file_stats.write("Number of genes with no chromosome info: " + str(gene_no_scaf_info)+ "\n")

out_file_stats.write("\nNumber of orths: " + str(N_orths ) + "\n")
out_file_stats.write("Number of orthologs Autosomal (soft): " + str(genes_A_soft_orth) + "\n")
out_file_stats.write("Number of orthologs  X-linked (soft): " + str(genes_X_soft_orth) + "\n")
out_file_stats.write("Number of orthologs  Autosomal (hard): " + str(genes_A_hard_orth) + "\n")
out_file_stats.write("Number of orthologs  X-linked (hard): " + str(genes_X_hard_orth) + "\n")



print("\n\nFinished, Walker.\n\n")	
 		