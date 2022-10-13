#!/usr/bin/python

#PYTHON SCRIPT 
#written by: Richard Wolfe
#
#to run type: python pullseq_header_name.py -i <inputfile> -o <outputfile> -n <file with names> -e <exclude T or F>
#         or: ./pullseq.py -i <inputfile> -o <outputfile> -m <min sequence length>
#
#   if error: /usr/bin/python^M: bad interpreter: No such file or directory
#      -there is a windows endl after shebang
#      -open in vi 
#         once in vi type:
#           :set ff=unix<return>
#           :x<return>
#
#
#
# extracts sequences from a fasta file 
#  matches to the first whitespace in the header line
#
#   -i <fasta file> (required) file to extract sequences from
#   -o <file name>  (required) name of file to write sequences into
#   -n <file name>  (required) file of header id names to search for
#   -e <T or F>     (required) exclude header names, T = exclude names, F = include names

import sys      #for exit command and maxint
import argparse #to get command line args 
                #needed to install argparse module because using python 2.6
                #and argparse comes with python 2.7
                #  sudo easy_install argparse


#create an argument parser object
#description will be printed when help is used
parser = argparse.ArgumentParser(description='A script to extract sequences from a fasta file')

#add the available arguments -h and --help are aqdded by default
#if the input file does not exist then program will exit
#if output file does not exit it will be created
# args.input is the input file Note: cant write to this file because read only
# args.output is the output file
# args.m is the minimum seq length
#Use fileType rU = universal so will read mac and windows endlines
parser.add_argument('-i', '--input', type=argparse.FileType('rU'), help='Input file name',required=True)
parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='Output file name', required=True)
parser.add_argument('-n', '--names', type=argparse.FileType('rU'), help='File of header id names', required=True)
parser.add_argument('-e', '--excluded', help='Exclude T or F', required=True)

#get the args
args = parser.parse_args()

#additional argument tests
if args.excluded != "T":
	if args.excluded != "F":
		print "Error: argument -e must be T or F"
		sys.exit(0)

#Test print the args
#print args

input_lines = 0
input_sequences = 0
output_lines = 0
output_sequences = 0


#read the ids into a list
#ids = args.names.readlines()

#remove the \n from each line in ids
#for index in xrange(len(ids)):
#	ids[index] = ids[index].rstrip() #removes whitespace from right of string
#close the file
#args.names.close()	


#read the ids to pull into a list
ids = []
ids_found = [] #keep track if id was found

line = args.names.readline()
while line:
	#print line
	line = line.strip() #removes whitespace from right of string
	line = line.split()[0] #only take 1st word in line 
	line = ">" + line
	ids.append(line)
	ids_found.append(0)
	line = args.names.readline()

args.names.close()



###############################
#read first line should be header line starts with ">"
line = args.input.readline()
input_lines = input_lines + 1

#if the file is not empty keep reading one at a time
while line:
	sequence = ""
	header = ""	

	#input_lines = input_lines + 1
	line = line.rstrip()  #remove whitespace at end includung endl
	header = line
	input_sequences = input_sequences + 1

	#read sequence
	line = args.input.readline()
	input_lines = input_lines + 1
	line = line.rstrip()  #remove whitespace at end includung endl
	

	while line:
		#input_lines = input_lines + 1
		#if line starts with > it is the header line for next sequence
		if line.startswith('>'):   #">" in line:
			#input_sequences = input_sequences + 1
			break   #break out of this while loop
				
		else: #this is the seq
			line = line.rstrip()  #remove whitespace at end includung endl
			sequence = sequence + line 
			#read another line
			line = args.input.readline()
			input_lines = input_lines + 1

	#the seq is found go through the ids and see ifthe header starts with an id
	found = False;
	header_split = header.split() #split the header on whitespace
	#for index in xrange(len(ids)):
	#	temp = ">" + ids[index]
        #        #print temp
	#	if header_split[0] == temp: #first char should be >
	#		found = True
	#		break
	
	if header_split[0] in ids:
		found = True
		#ids.remove(header_split[0]) #if duplicates then will only select 1
		index = ids.index(header_split[0])
		ids_found[index] += 1  #increment the found count

	if found:
		if args.excluded == "F":
			args.output.write(header + "\n")
			args.output.write(sequence + "\n")
			output_sequences = output_sequences + 1
	else:
		if args.excluded == "T":
 			args.output.write(header + "\n")
			args.output.write(sequence + "\n")
			output_sequences = output_sequences + 1

		#else:  #should not get here
			#print "ERROR......args.excluded not T or F"
			#sys.exit(1)



#close the files
args.input.close()
args.output.close()


print "Lines read from input file = ", input_lines
print "Sequences in input file = ", input_sequences
print "Sequences written to output file = ", output_sequences


#print if found data
no_ids_found = 0
more_than_one_found = 0


for item in ids_found:
	if item == 0:
		no_ids_found += 1
	elif item > 1:
		more_than_one_found += 1

print ""
print "The number of sequences not found in the fasta file = ", no_ids_found
print "The number of sequences written to the outputfile that have the same header = ",more_than_one_found

print ""
print "Script finished..."
