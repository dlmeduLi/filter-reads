#!/usr/bin/python

from __future__ import print_function

import re
import sys
import os
import os.path
from optparse import OptionParser
import pysam
import csv

# sam regular expressions

mdRe = re.compile('^[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*$')
mdTagRe = re.compile('[0-9]+|\^[A-Z]+|[A-Z]')
pairKeyRe = re.compile('(.*)\s[1-2](.*)')

# Count file lines

def opcount(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

# Get alignment unmatched bases

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

	return um

def main():
	# parse the command line options
	
	usage = 'usage: %prog [options] input.sam snp1.csv snp2.csv'
	parser = OptionParser(usage=usage, version='%prog version 0.0.1')

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

	samBaseName = os.path.splitext(os.path.basename(samFileName))[0]
	outputSamFile = samBaseName + '.undetermined.sam'
	outputSamFile1 = samBaseName + '.' + os.path.splitext(os.path.basename(snpFileName1))[0] + '.sam'
	outputSamFile2 = samBaseName + '.' + os.path.splitext(os.path.basename(snpFileName2))[0] + '.sam'

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

	# build pair tables

	pairSnp1 = {}
	pairSnp2 = {}

	print('* Processing...')
	for read in samfile.fetch():

		pairKey = ''.join(pairKeyRe.findall(read.qname)[0])

		# look up pair tables

		if(pairKey in pairSnp1):
			outSam1.write(read)
			writtenLineCount1 += 1
			pairSnp1.pop(pairKey, None)
			continue

		if(pairKey in pairSnp2):
			outSam2.write(read)
			writtenLineCount2 += 1
			pairSnp2.pop(pairKey, None)
			continue
			
		# split

		if (not read.has_tag('MD')) :
			outSam.write(read)
			writtenLineCount += 1
			continue

		unmatches = GetUnmatchList(read.seq, read.cigar, read.get_tag('MD'))
		chrname = samfile.getrname(read.rname)
		basepos = read.pos + 1 # the returned pos is 0 based

		# select alignments 

		issnp1 = False
		issnp2 = False
		for um in unmatches :
			key = chrname + ':' + str(um[0] + basepos) + ':' + um[1]

			# look up snp tables

			issnp1 = ((key in snp1) and (snp1[key] == um[2]))
			issnp2 = ((key in snp2) and (snp2[key] == um[2]))

			if(issnp1 and (not issnp2)):
				outSam1.write(read)
				writtenLineCount1 += 1

				# store pair key

				if(read.is_paired):
					pairSnp1[pairKey] = True

				break
			elif(issnp2 and (not issnp1)):
				outSam2.write(read)
				writtenLineCount2 += 1
				
				# store pair key

				if(read.is_paired):
					pairSnp2[pairKey] = True

				break
		
		if((not issnp1) and (not issnp2)):
			outSam.write(read)
			writtenLineCount += 1

		# progress 

		lineCount = lineCount + 1
		if(totalLineCount == 0):
			percentage = 0
		else:
			percentage = lineCount * 1.0 / totalLineCount
		sys.stdout.write('\r  read %ld (%.2f%%), (snp1: %ld, snp2: %ld, undt: %ld)' 
						% (lineCount, percentage * 100, 
						writtenLineCount1, writtenLineCount2, writtenLineCount))
		sys.stdout.flush()

	sys.stdout.write('\r  read %ld (%.2f%%)' % (lineCount, 100))
	sys.stdout.flush()
	
	samfile.close()
	
	# Clear resources

	outSam.close()
	outSam1.close()
	outSam2.close()

	print('\n* Complete')
	readCount = writtenLineCount + writtenLineCount1 + writtenLineCount2
	if(readCount == 0):
		snp1Percent = 0.0
		snp2Percent = 0.0
		undtPercent = 0.0
	else:
		snp1Percent = writtenLineCount1 / readCount
		snp2Percent = writtenLineCount2 / readCount
		undtPercent = writtenLineCount / readCount
	print('  %ld reads (%.2f%%) written to %s.\n  %ld reads (%.2f%%) written to %s.\n  %ld reads (%.2f%%) written to %s.' 
		  % (writtenLineCount, undtPercent * 100, os.path.basename(outputSamFile),
		  	 writtenLineCount1, snp1Percent * 100, os.path.basename(outputSamFile1),
		  	 writtenLineCount2, snp2Percent * 100, os.path.basename(outputSamFile2)))


if __name__ == '__main__':
	main()