#!/usr/bin/python

from __future__ import print_function

import re
import sys
import os
import os.path
import subprocess
from optparse import OptionParser
import shelve
import pysam
import csv

# sam regular expressions

mdRe = re.compile('^[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*$')
mdTagRe = re.compile('[0-9]+|\^[A-Z]+|[A-Z]')

# Sort SAM file with system sort program
# sed -i '/^@/d;/!^@/q' file

def SortSamFile(inputFile, outputFile):
	if(not os.path.exists(inputFile)):
		return False

	# Extract the header of the input SAM file

	os.system('sed -n \'/^@/p;/^[^@]/q\' ' + inputFile + ' > ' + outputFile)

	# Remove the header in the temp file

	os.system('sed -n \'/^@/d;/^[^@]/p\' ' + inputFile + ' | sort >> ' + outputFile)
	
	return True

# Count file lines

def opcount(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

# Get the qname tag

def QNameKey(qname, keyRe):
	key = qname
	if keyRe :
		if(keyRe.match(qname)):
			tags = keyRe.findall(key)
			if(len(tags) > 0):
				key = ''
				if(type(tags[0]) is str):
					key = tags[0]
				elif(type(tags[0]) is tuple):
					for s in tags[0] :
						key += s
	return key

# Calculate the key an alignment
# The algrithm should make sure that one key identify one unique read

def AlignmentKey(alignment, keyRe):
	# Get the tag of QNAME without read number
	key = QNameKey(alignment.qname, keyRe)

	if(alignment.flag & 0x40):
		key += (':' + str(alignment.pos) + ':' + str(alignment.pnext))
	elif(alignment.flag & 0x80):
		key += (':' + str(alignment.pnext) + ':' + str(alignment.pos))

	return key

def AlignmentGroupKey(alignment, keyRe):
	return QNameKey(alignment.qname, keyRe)

def GetUnmatchList(seq, cigar, md):
	um = []
	refseq = ''
	pos = 0
	for m in cigar :
		if(m[0] == 0) : # M
			refseq += seq[pos : pos + m[1]]
			pos += m[1]
		elif(m[0] == 1): # I
			pos += m[1]
		elif(m[0] == 2): # D
			refseq += ('N' * m[1])
		elif(m[0] == 3): # N
			refseq += ('N' * m[1])
		elif(m[0] == 4): # S
			pos += m[1]

	# Parse MD & extract the unmatched bases

	pos = 0
	tags = mdTagRe.findall(md)
	try:
		for tag in tags:
			if(tag.isdigit()):
				pos += int(tag)
			elif(tag.isalpha()):

				# Get the base in the real reads

				um += [(pos, tag, refseq[pos])]
				pos += 1
			elif(tag[0] == '^'):
				pos += (len(tag) - 1)
			else:
				pos += len(tag)
	except IndexError:
		return []

	# print(seq)
	# print(refseq)
	return um

def main():
	# parse the command line options
	
	usage = 'usage: %prog [options] input.sam snp1.csv snp2.csv'
	parser = OptionParser(usage=usage, version='%prog version 0.0.1')
	
	# parser.add_option('-o', '--output-file', dest='outputfile',
	# 					help='write the result to output file')
	# parser.add_option('-s', '--sort', 
	# 					action="store_true", dest="sort", default=False,
	# 					help='sort the input SAM file before further processing')
	# parser.add_option('-k', '--key-reg', dest="keyreg",
	# 					help='qname regular expression to extract the alignment key')

	(options, args) = parser.parse_args()
	if(len(args) != 3):
		parser.print_help()
		sys.exit(0)

	samFileName = args[0]
	snpFileName1 = args[1]
	snpFileName2 = args[2] 

	if(not os.path.exists(samFileName)):
		print('error: Failed to open file "', samFileName, '"')
		sys.exit(-1)

	if(not os.path.exists(snpFileName1)):
		print('error: Failed to open file "', snpFileName1, '"')
		sys.exit(-1)

	if(not os.path.exists(snpFileName2)):
		print('error: Failed to open file "', snpFileName2, '"')
		sys.exit(-1)

	print('* Initializing...')

	outputSamFile = "uncertain.sam"
	outputSamFile1 = os.path.splitext(os.path.basename(snpFileName1))[0] + '.sam'
	outputSamFile2 = os.path.splitext(os.path.basename(snpFileName2))[0] + '.sam'

	samfile = pysam.AlignmentFile(samFileName, 'r')
	outSam = pysam.AlignmentFile(outputSamFile, 'wh', template = samfile, text = samfile.text)
	outSam1 = pysam.AlignmentFile(outputSamFile1, 'wh', template = samfile, text = samfile.text)
	outSam2 = pysam.AlignmentFile(outputSamFile2, 'wh', template = samfile, text = samfile.text)

	totalLineCount = opcount(samFileName)
	print('  %ld lines.' % totalLineCount)
	lineCount = 0
	writtenLineCount = 0
	writtenLineCount1 = 0
	writtenLineCount2 = 0

	# load SNP data sets

	print('  loading SNP data file \"%s\" ...' % (os.path.basename(snpFileName1)))
	snp1 = {}
	csvSnp1 = csv.reader(open(snpFileName1))
	for row in csvSnp1 :
		snp1['chr' + str(row[0]) + ':' + str(row[1]) + ':' + str(row[2])] = row[3]																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																
	
	print('  loading SNP data file \"%s\" ...' % (os.path.basename(snpFileName2)))
	snp2 = {}
	csvSnp2 = csv.reader(open(snpFileName2))
	for row in csvSnp2 :
		snp2['chr' + str(row[0]) + ':' + str(row[1]) + ':' + str(row[2])] = row[3]

	print('* Processing...')
	for read in samfile.fetch():

		if (not read.has_tag) :
			outSam.write(read)
			writtenLineCount += 1
			continue

		unmatches = GetUnmatchList(read.seq, read.cigar, read.get_tag('MD'))
		chrname = samfile.getrname(read.rname)
		basepos = read.pos

		# select alignments 

		found = False
		for um in unmatches :
			key = chrname + ':' + str(um[0] + basepos) + ':' + um[1]

			# look up table snp1

			if((key in snp1) and (snp1[key] == um[2])):
				outSam1.write(read)
				writtenLineCount1 += 1
				found = True
				break
			elif((key in snp2) and (snp2[key] == um[2])):
				outSam2.write(read)
				writtenLineCount2 += 1
				found = True
				break 
		
		if(not found):
			outSam.write(read)
			writtenLineCount += 1

		# progress 

		lineCount = lineCount + 1
		if(totalLineCount == 0):
			percentage = 0
		else:
			percentage = lineCount * 1.0 / totalLineCount
		sys.stdout.write('\r  read %ld (%.2f%%)' % (lineCount, percentage * 100))
		sys.stdout.flush()

	sys.stdout.write('\r  read %ld (%.2f%%)' % (lineCount, 100))
	sys.stdout.flush()
	
	samfile.close()
	
	# Clear resources

	outSam.close()
	outSam1.close()
	outSam2.close()

	print('\n* Complete')
	print('  %ld reads written to %s.\n  %ld reads written to %s.\n  %ld reads written to %s.' 
		  % (writtenLineCount, os.path.basename(outputSamFile),
		  	 writtenLineCount1, os.path.basename(outputSamFile1),
		  	 writtenLineCount2, os.path.basename(outputSamFile2)))


if __name__ == '__main__':
	main()