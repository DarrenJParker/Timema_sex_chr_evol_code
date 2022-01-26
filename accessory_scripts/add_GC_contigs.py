import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'c:g:o:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

GC_file_name  = None
chr_file_name = None
outprefix     = "testout"

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** default_py.py | Written by DJP, 02/04/20 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Template script")	
		print("\n**** Usage****\n")
		print("python3 default_py.py -i [in file] -o [outprefix] \n\n")
		sys.exit(2)
	elif opt in ('-g'):
		GC_file_name  = arg
	elif opt in ('-c'):
		chr_file_name  = arg
	elif opt in ('-o'):
		outprefix = arg
	else:
		print("i dont know")
		sys.exit(2)


### read in GC


GC_dict = {}
GC_file = open(GC_file_name)
for line in GC_file:
	line = line.rstrip("\n").split(",")
	GC_dict[line[0]] = line[1]

GC_file.close()


### add to chr info

outfile = open(chr_file_name.replace(".csv", "_wGC.csv"), "w")

chr_file = open(chr_file_name)
line_N = 0
for line in chr_file:
	line_N = line_N + 1
	line = line.rstrip("\n").replace('"', '')
	#print(line)
	if line_N == 1:
		outfile.write(line + ",GC\n" )
	else:
		scaf = line.rstrip("\n").split(",")[0]
		GC_v = GC_dict.get(scaf)
		outfile.write(line + "," + GC_v + "\n" )

chr_file.close()



