from __future__ import print_function
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(prog='tbdream.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-gff', '--gff', type=str, help='Gff file with genes annotation', required=True)
parser.add_argument('vcf', type=str, help='Vcf file with snps')
#parser.add_argument('-td', '--trainingDockedMol2', metavar='GlobConfig', type=str, help='Full path to docked training multiMol2', required=True)
#parser.add_argument('-bd', '--baseDockedMol2', metavar='GlobConfig', type=str, help='Full path to database docked multiMol2', required=True)


args = parser.parse_args()

SUPRESSWARNINGS = False

class geneRecord():
	start = 0
	end = 0
	strand = 1
	description = ''
	gid = ''
	def __init__(self, gid, start, end, strand, description):
		self.gid = gid
		self.start = start
		self.end = end
		self.strand = strand
		self.description = description

def processGffGenes(fileName):
	genes = []
	fHandle = open(fileName)
	Crtype = 2
	Cgstart = 3
	Cgend = 4
	Cgstrand = 6
	Cgdescription = 8
	for line in fHandle:
		arr = line.strip().split()
		if len(arr):
			gid = ""
			for i in arr[Cgdescription].strip().split(';'):
				if "ID=" in i:
					gid = i.strip()[3:].replace("RVBD_", 'Rv')
			if arr[Cgstrand].strip() == '+':
				strand = 1
			else:
				strand = -1
			genes.append(geneRecord(gid, int(arr[Cgstart].strip()), int(arr[Cgend].strip()), strand, arr[Cgdescription]))
	fHandle.close()
	return genes

def getGene(genes, geneId):
	for gene in genes:
		if geneId in gene.gid:
			return gene
	return None
class SNV():
	stype = ''
	position = 0
	ref = ''
	alt = ''
	gid = ''
	pubmed = ''
	drug = ''
	bid = ''
	def __init__(self, stype, position, ref, alt, gid, pubmed, drug, bid = ''):
		self.stype = stype
		self.position = position
		self.ref = ref
		self.alt = alt
		self.gid = gid
		self.pubmed = pubmed
		self.drug = drug
		self.bid = bid


reverce={}
reverce['A'] = 'T'
reverce['T'] = 'A'
reverce['G'] = 'C'
reverce['C'] = 'G'

def processDatabase1(fileName, genes):
	sorted(genes, key=lambda x: x.start)
	snvs = {}
	fHandle = open(fileName)
	Cgid = 8
	Cdrug = 9
	Cref = 7
	Cnuclpos = 1
	Ctype = 2
	Crefn = 3
	Caltn = 4
	Caarepl = 7
	Cpubmed = 10
	header = fHandle.next()
	for line in fHandle:
		items = line.strip().split('\t')
		if len(items) < 9:
			continue
		geneId = items[Cgid]
		gene = getGene(genes, geneId)
		if not gene:
			if not SUPRESSWARNINGS:
				print("There is no gene " + geneId)
			continue
		if gene.strand > 0:
			genomePosition = gene.start + int(items[Cnuclpos]) - 1
		else:
			genomePosition = gene.end - int(items[Cnuclpos])
		refS = ""
		altS = ""
		for i in items[Crefn]:
			if i in "ATGC":
				if gene.strand > 0:
					refS += i
				else:
					refS += reverce[i]
		for i in items[Caltn]:
			if i in "ATGC":
				if gene.strand > 0:
					altS += i
				else:
					altS += reverce[i]
		entryType = items[Ctype]
		entryName = str(genomePosition) + '_'
		if entryType == 'del':
			entryName += 'del_' + refS
		elif entryType == 'ins':
			entryName += 'ins_' + altS
		else:
			entryName += refS + '_' + altS
		snvs[entryName] =  SNV(items[Ctype], genomePosition, refS, altS, geneId, items[Cpubmed], items[Cdrug], items[0])
	fHandle.close()
	return snvs

class geneInfo():
	gid = ''
	descr = ''
	start = 0
	end = 0
	strand = 1
	seq = ''
	def __init__(self, gid, start, end, strand, seq = '', descr = '', ):
		self.gid = gid
		self.start = start
		self.end = end
		self.strand = strand
		self.seq = seq
		self.descr = descr


class tbvarSNV():
	stype = ''
	position = 0
	ref = ''
	alt = ''
	gene = ''
	freq = 0
	siftPred = ''
	def __init__(self, stype, position, ref, alt, gene, freq, siftPred):
		self.stype = stype
		self.position = position
		self.ref = ref
		self.alt = alt
		self.gene = gene
		self.freq = freq
		self.siftPred = siftPred


def processDatabase2(fileName, genes):
	sorted(genes, key=lambda x: x.start)
	snvs = {}
	fHandle = open(fileName)
	Cnuclpos = 0
	Crefn = 2
	Caltn = 3
	Ctype = 5
	Cfreq = 6
	Cgid = 7
	CgeneName = 8
	CgeneStart = 9
	CgeneEnd = 10
	CgeneStrand = 11
	CsiftPred = 13

	header = fHandle.next()
	for line in fHandle:
		items = line.strip().split('\t')
		if CsiftPred > len(items) or not len(items[Cgid]):
			continue
		geneId = items[Cgid]
		gene = geneInfo(geneId, items[CgeneStart], items[CgeneEnd], items[CgeneStrand])
		genomePosition = int(items[Cnuclpos])

		refS = items[Crefn].strip().replace(' ', '')
		altS = items[Caltn].strip().replace(' ', '')

		snvs[str(genomePosition) + '_' + refS + '_' + altS] = tbvarSNV(items[Ctype], genomePosition, refS, altS, gene, items[Cfreq], items[CsiftPred])
	fHandle.close()
	return snvs

def getReferenceFasta(fastaFileName):
	handle = open(fastaFileName)
	fasta = None
	for record in SeqIO.parse(handle, "fasta"):
		fasta=record
		break
	handle.close()
	return fasta

class SNP():
	chrom = ''
	ref = ''
	alt = ''
	pos = 0
	def __init__(self, chrom, pos, ref, alt):
		self.chrom = chrom
		self.pos = pos
		self.ref = ref
		self.alt = alt

def readVcf(fileName):
	handle = open(fileName)
	snps = {}
	Cchrom = 0
	Cpos = 1
	Cref = 3
	Calt = 4
	#rest is not informative for now
	for line in handle:
		if line.startswith('##'):
			pass
		elif line.startswith('#'):
			header = line.strip().split()
		else:
			items = line.strip().split()
			pos = int(items[Cpos])
			entryName = str(pos) + '_' + items[Cref] + '_' + items[Calt]
			snps[entryName] = SNP(items[Cchrom], pos, items[Cref], items[Calt])
	handle.close()
	return snps
def searchSnps(snvs, vcfSnps):
	sorted(snvs, key=lambda x: x.position)
	res = {}
	for snpPos in vcfSnps:
		filtered = []
		for snv in snvs:
			if snv.position == snpPos:
				filtered.append(snv)
		res[snpPos] = filtered
	return res


def getGeneSeq(fasta, gene):
	print(len(fasta), gene.start, gene.end)
	subseq = fasta[gene.start: gene.end]
	if gene.strand < 0:
		subseq = subseq.reverse_complement()
	return subseq
genes = processGffGenes(args.gff)
referenceFasta = getReferenceFasta("H37RV_V5.fasta")#FIXME
dreamSnvs = processDatabase1("dbTest", genes)
tbvarSnvs = processDatabase2("base2Annotation", genes)

vcfSnps = readVcf(args.vcf)
for snp in vcfSnps:
	print(snp, end = '', sep = '\t')
	if snp in dreamSnvs:
		print('\t!', dreamSnvs[snp].gid, dreamSnvs[snp].drug, dreamSnvs[snp].pubmed, sep = '\t')
	elif snp in tbvarSnvs:
		print('\t?', tbvarSnvs[snp].gene.gid, tbvarSnvs[snp].freq, tbvarSnvs[snp].siftPred, sep = '\t')
	else:
		print('\t.\tno_match')

# for snp in dreamSnvs:
	# print(snp, dreamSnvs[snp].bid)

#foundDream = searchSnps(snvs, dreamSnvs, tbvarSnvs)
# for snp in vcfSnps:
# 	print(vcfSnps[snp].pos, vcfSnps[snp].ref, vcfSnps[snp].alt, len(found[snp]))
# 	for i in found[snp]:
# 		print('Dream', i.drug, i.ref, i.alt)

# for i in snvs:
# 	print(i.position, i.stype, i.ref, i.alt, i.geneId, i.pubmed, i.drug)
# for gene in genes:
# 	if "1908c" in gene.gid:
# 		s = getGeneSeq(referenceFasta, gene).seq
# 		break
# print(s[940:945])
#sorted(genes) TODO
#for i in genes:
#	print(i.gid, i.start, i.end)
