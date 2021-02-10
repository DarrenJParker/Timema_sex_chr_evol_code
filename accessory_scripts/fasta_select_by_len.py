## fasta_select_by_len 

import sys
import os

args = sys.argv
arg_len = len(args)

if arg_len <5:
	print("\n**** Written by DJP, 24/11/15 in Python 3.4 ****\n")
	print("This program takes a fasta file as input.") 
	print("It outputs a fasta file filtered by specified length.")
	print("\n**** USAGE **** \n")
	print("fasta_select_by_len.py [name of fasta file] [max len] [min len] [output file]\n")
	print("\n**** EXAMPLE ****\n")
	

	print("to output sequences greater-or-equal-to than 200bp:")
	print("fasta_select_by_len.py [name of fasta file] [max] [200] [output file]\n")
	
	print("to output sequences less-or-equal-to than 200bp:")
	print("fasta_select_by_len.py [name of fasta file] [200] [min] [output file]\n")

	print("to output sequences between (and including) 400 and 200bp:")
	print("fasta_select_by_len.py [name of fasta file] [400] [200] [output file]\n")



else:
	input_fasta = args[1]
	greater_than = args[2]
	less_than = args[3]
	out_file = args[4]
	
	##### FIRST unwrap fasta - precautionary will be necessary for some files 
	### note making a temp unwrapped fasta file  - removed at end

	output_fasta_name = input_fasta + ".TEMP_extract_fasta_file" 

	output_file = open(output_fasta_name, "w")
	print("\nUnwrapping fasta file")
	count = 0
	in_file = open(input_fasta)
	for line in in_file:
		count = count + 1
		line = line.rstrip("\n")
		if line.startswith(">") and count == 1:
			output_file.write(line + "\n")
		elif line.startswith(">") and count > 1:
			output_file.write("\n" + line + "\n")
		else: 
			output_file.write(line)	

	output_file.close()


	### add seqs to dictionary
	name_list = []
	seq_list = []
	seq_dict = {}
	
	done = 0
	seq_file_1 = open(output_fasta_name)
	for line in seq_file_1:
		lineA = line.rstrip("\n")
		if lineA.startswith(">"):
			lineB = lineA.replace(">", "")
			name_list.append(lineB)
		else:
			seq_list.append(lineA)
			done = done + 1
			# done_divide = done / 1000
			# if done_divide.is_integer():
			# 	print("Read " + str(done) + " sequences from " + input_fasta)

	for element in range(0,len(name_list)):
		name1 = name_list[element]
		seq1 = seq_list[element]
		seq_dict[name1] = seq1
		
	print("Read " + str(len(name_list)) + " sequences from " + input_fasta)

	#print(seq_dict)
	## tidyup
	seq_file_1.close()
	os.remove(output_fasta_name)
	
	output_N_seq = 0
	
	output_file = open(out_file,"w")
	if less_than != "min" and greater_than == "max":
		print("\nOutputting sequences greater-than-or-equal-to " + str(less_than) + "\n")
		for el in name_list:
			seq_len = len(seq_dict.get(el))
			#print(seq_len)
			if seq_len >= int(less_than):
				output_file.write(">" + el + "\n")
				output_file.write(seq_dict.get(el) + "\n")
				output_N_seq = output_N_seq + 1
						
	elif less_than == "min" and greater_than != "max":
		print("\nOutputting sequences less-than-or-equal-to " + str(greater_than) + "\n")
		for el in name_list:
			seq_len = len(seq_dict.get(el))
			#print(seq_len)
			if seq_len <= int(greater_than):
				output_file.write(">" + el + "\n")
				output_file.write(seq_dict.get(el) + "\n")
				output_N_seq = output_N_seq + 1

	elif less_than != "min" and greater_than != "max":
		print("\nOutputting sequences between (and including) " + str(greater_than) + " and " + str(less_than) + "\n")
		for el in name_list:
			seq_len = len(seq_dict.get(el))
			#print(seq_len)
			if seq_len <= int(greater_than) and seq_len >= int(less_than):
				output_file.write(">" + el + "\n")
				output_file.write(seq_dict.get(el) + "\n")
				output_N_seq = output_N_seq + 1
				
	else:
		print("\nPlease specify max and min values!\n")

	
	print("Total seqs outputted to " + out_file +  ": " + str(output_N_seq) + "\n\n")
