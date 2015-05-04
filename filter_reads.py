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
tagCTRe = re.compile('.*CT.*')
tagGARe = re.compile('.*GA.*')

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

def IsUnmatchesInSnps(unmatches, snps, chrname, basepos):
	for um in unmatches:
		key = chrname + ':' + str(um[0] + basepos)
		value =  um[1] + ':' + um[2]
		if((key in snps) and (value == snps[key])):
			return True
	return False

def HandleReadUndt(key, read, samfile, dictUndt):
	# check if there is paired read

	if(key in dictUndt):
		samfile.write(dictUndt[key])
		dictUndt.pop(key, None)
		samfile.write(read)

		return 2

	# unpaired read, write directly to undetermined file

	if(not read.is_proper_pair):
		samfile.write(read)

		return 1

	# paired read, save it for future lookup

	dictUndt[key] = read

	return 0

def HandleReadSnp(key, read, samfile, dictSnp, dictUndt):
	# check if there is paired read

	if(key in dictUndt):
		samfile.write(dictUndt[key])
		dictUndt.pop(key, None)
		samfile.write(read)

		return 2

	# unpaired read, write directly to snp1 file

	if(not read.is_proper_pair):
		samfile.write(read)

		return 1

	# paired read, write read & save the key for future lookup

	samfile.write(read)
	dictSnp[key] = True

	return 1

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
	samFileNameUndt = samBaseName + '.undetermined.sam'
	samFileNameSnp1 = samBaseName + '.' + os.path.splitext(os.path.basename(snpFileName1))[0] + '.sam'
	samFileNameSnp2 = samBaseName + '.' + os.path.splitext(os.path.basename(snpFileName2))[0] + '.sam'
	logFileName = samBaseName + '.log'

	samfile = pysam.AlignmentFile(samFileName, 'r')
	samFileUndt = pysam.AlignmentFile(samFileNameUndt, 'wh', template = samfile, text = samfile.text)
	samFileSnp1 = pysam.AlignmentFile(samFileNameSnp1, 'wh', template = samfile, text = samfile.text)
	samFileSnp2 = pysam.AlignmentFile(samFileNameSnp2, 'wh', template = samfile, text = samfile.text)
	try:
		logfile = open(logFileName, 'w')
	except IOError:
		print('error: create log file failed!')
		sys.exit(-1)

	totalLineCount = opcount(samFileName)
	print('  %ld lines.' % totalLineCount)
	lineCount = 0
	countUndt = 0
	countSnp1 = 0
	countSnp2 = 0

	# load SNP data sets

	print('  loading SNP data file \"%s\" ...' % (os.path.basename(snpFileName1)))
	dataSnp1 = {}
	csvSnp1 = csv.reader(open(snpFileName1))
	for row in csvSnp1 :
		dataSnp1['chr' + str(row[0]) + ':' + str(row[1])] = row[2] + ':' + row[3]																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																
	
	print('  loading SNP data file \"%s\" ...' % (os.path.basename(snpFileName2)))
	dataSnp2 = {} 
	csvSnp2 = csv.reader(open(snpFileName2))
	for row in csvSnp2 :
		dataSnp2['chr' + str(row[0]) + ':' + str(row[1])] = row[2] + ':' + row[3]

	# build pair tables

	dictUndt = {}
	dictSnp1 = {}
	dictSnp2 = {}

	print('* Processing...')
	for read in samfile.fetch():

		chrname = samfile.getrname(read.rname)
		pairKey = ''.join(pairKeyRe.findall(read.qname)[0])
		tagXB = read.get_tag('XB')
		isCT = tagCTRe.match(tagXB)
		isGA = tagGARe.match(tagXB)

		# look up pair tables

		if(pairKey in dictSnp1):
			samFileSnp1.write(read)
			countSnp1 += 1
			dictSnp1.pop(pairKey, None)

			lineCount += 1
			continue

		if(pairKey in dictSnp2):
			samFileSnp2.write(read)
			countSnp2 += 1
			dictSnp2.pop(pairKey, None)

			lineCount += 1
			continue

		# Get the read covered SNPS

		dictCoveredSnp1 = {}
		dictCoveredSnp2 = {}
		refPoses = read.get_reference_positions()
		for pos in refPoses:
			posKey = chrname + ':' + str(pos + 1)
			if(posKey in dataSnp1):
				snp = dataSnp1[posKey]
				if isCT and (snp == 'C:T' or snp == 'T:C'):
					continue
				elif isGA and (snp == 'G:A' or snp == 'A:G'):
					continue
				dictCoveredSnp1[posKey] = dataSnp1[posKey]
			if(posKey in dataSnp2):
				if isCT and (snp == 'C:T' or snp == 'T:C'):
					continue
				elif isGA and (snp == 'G:A' or snp == 'A:G'):
					continue
				dictCoveredSnp2[posKey] = dataSnp2[posKey]

		isCoveredSnp1 = (len(dictCoveredSnp1) > 0)
		isCoveredSnp2 = (len(dictCoveredSnp2) > 0)

		# Get Mismatches

		unmatches = GetUnmatchList(read.seq, read.cigar, read.get_tag('MD'))
		hasUnmatches = (len(unmatches) > 0)


		basepos = read.pos + 1 # the returned pos is 0 based
		isUnmatchesInSnp1 = IsUnmatchesInSnps(unmatches, dictCoveredSnp1, chrname, basepos)
		isUnmatchesInSnp2 = IsUnmatchesInSnps(unmatches, dictCoveredSnp2, chrname, basepos)

		# Split Strategy

		if((not isCoveredSnp1) and (not isCoveredSnp2)):
			
			# Undetermined

		 	countUndt += HandleReadUndt(pairKey, read, samFileUndt, dictUndt)

		elif((isCoveredSnp1) and (not isCoveredSnp2)):
			if(not hasUnmatches):
				
				# Snp2

				countSnp2 += HandleReadSnp(pairKey, read, samFileSnp2, dictSnp2, dictUndt)

			else:
				if(isUnmatchesInSnp1):

					# Snp1

					countSnp1 += HandleReadSnp(pairKey, read, samFileSnp1, dictSnp1, dictUndt)

				else:

					# Undetermined

					countUndt += HandleReadUndt(pairKey, read, samFileUndt, dictUndt)

		elif((not isCoveredSnp1) and isCoveredSnp2):
			if(not hasUnmatches):

				# Snp1

				countSnp1 += HandleReadSnp(pairKey, read, samFileSnp1, dictSnp1, dictUndt)
			else:
				if(isUnmatchesInSnp2):

					# Snp2

					countSnp2 += HandleReadSnp(pairKey, read, samFileSnp2, dictSnp2, dictUndt)
				else:

					# Undetermined

					countUndt += HandleReadUndt(pairKey, read, samFileUndt, dictUndt)

		elif(isCoveredSnp1 and isCoveredSnp2):
			if(isUnmatchesInSnp1 and (not isUnmatchesInSnp2)):

				# Snp1

				countSnp1 += HandleReadSnp(pairKey, read, samFileSnp1, dictSnp1, dictUndt)

			elif((not isUnmatchesInSnp1) and isUnmatchesInSnp2):

				# Snp2

				countSnp2 += HandleReadSnp(pairKey, read, samFileSnp2, dictSnp2, dictUndt)

			else:

				# Undetermined

				countUndt += HandleReadUndt(pairKey, read, samFileUndt, dictUndt)
		else:
			# Undetermined

		 	countUndt += HandleReadUndt(pairKey, read, samFileUndt, dictUndt)

		# progress 

		lineCount = lineCount + 1
		if(totalLineCount == 0):
			percentage = 0
		else:
			percentage = lineCount * 1.0 / totalLineCount
		sys.stdout.write('\r  read #%ld (%.2f%%), (SNP1: %ld, SNP2: %ld, UNDT: %ld)' 
						% (lineCount, percentage * 100, 
						countSnp1, countSnp2, countUndt))
		sys.stdout.flush()

	# Write back the undetermined dict records

	for key, read in dictUndt.iteritems() :
		countUndt += 1
		samFileUndt.write(read)

	logfile.write('\nUnpaired read count: %ld\n' %(len(dictUndt)))

	sys.stdout.write('\r  read #%ld (%.2f%%)' % (lineCount, 100))
	sys.stdout.flush()
	
	# Clear resources

	samfile.close()
	samFileUndt.close()
	samFileSnp1.close()
	samFileSnp2.close()

	print('\n* Complete')
	countTotal = countUndt + countSnp1 + countSnp2
	if(countTotal == 0):
		percentSnp1 = 0.0
		percentSnp2 = 0.0
		percentUndt = 0.0
	else:
		percentSnp1 = countSnp1 * 1.0 / countTotal
		percentSnp2 = countSnp2 * 1.0 / countTotal
		percentUndt = countUndt * 1.0 / countTotal
	stats = '  %ld reads (%.2f%%) written to %s.\n  %ld reads (%.2f%%) written to %s.\n  %ld reads (%.2f%%) written to %s. \n  total: %ld reads.' % (
			countUndt, percentUndt * 100.0, os.path.basename(samFileNameUndt),
		  	countSnp1, percentSnp1 * 100.0, os.path.basename(samFileNameSnp1),
		  	countSnp2, percentSnp2 * 100.0, os.path.basename(samFileNameSnp2),
		  	countTotal)
	print(stats)
	logfile.write(stats)
	logfile.close()

if __name__ == '__main__':
	main()