import argparse
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(prog='tbdream.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-gff', '--gff', type=str, help='Full path to database multiMol2', required=True)
#parser.add_argument('-td', '--trainingDockedMol2', metavar='GlobConfig', type=str, help='Full path to docked training multiMol2', required=True)
#parser.add_argument('-bd', '--baseDockedMol2', metavar='GlobConfig', type=str, help='Full path to database docked multiMol2', required=True)


args = parser.parse_args()

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
					gid = i.strip()[3:]
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
	geneId = ''
	pubmed = ''
	drug = ''
	def __init__(self, stype, position, ref, alt, geneId, pubmed, drug):
		self.stype = stype
		self.position = position
		self.ref = ref
		self.alt = alt
		self.geneId = geneId
		self.pubmed = pubmed
		self.drug = drug

reverce={}
reverce['A'] = 'T'
reverce['T'] = 'A'
reverce['G'] = 'C'
reverce['C'] = 'G'

def processDatabase(fileName, genes):
	sorted(genes, key=lambda x: x.start)
	snvs = []
	fHandle = open(fileName)
	Cgid = 1
	Cdrug = 3
	Cref = 7
	Cnuclpos = 10
	Ctype = 11
	Crefn = 12
	Caltn = 13
	Caarepl = 14
	Cpubmed = 30
	header = fHandle.next()
	for line in fHandle:
		items = line.strip().split('\t')
		if not len(items):
			continue
		geneId = 'RVBD_' + items[Cgid]
		gene = getGene(genes, geneId)
		if not gene:
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
		snvs.append(SNV(items[Ctype], genomePosition,refS, altS, geneId, items[Cpubmed], items[Cdrug]))
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

def getGeneSeq(fasta, gene):
	print(len(fasta), gene.start, gene.end)
	subseq = fasta[gene.start: gene.end]
	if gene.strand < 0:
		subseq = subseq.reverse_complement()
	return subseq
genes = processGffGenes(args.gff)
referenceFasta = getReferenceFasta("H37RV_V5.fasta")#FIXME
snvs = processDatabase("dbTest", genes)
sorted(snvs, key=lambda x: x.position)
for i in snvs:
	print(i.position - 1, i.stype, i.ref, i.alt, i.geneId, i.pubmed, i.drug)
# for gene in genes:
# 	if "1908c" in gene.gid:
# 		s = getGeneSeq(referenceFasta, gene).seq
# 		break
# print(s[940:945])
#sorted(genes) TODO
#for i in genes:
#	print(i.gid, i.start, i.end)
