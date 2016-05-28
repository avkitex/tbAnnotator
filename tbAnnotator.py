from __future__ import print_function
import argparse, os, json
#from Bio import SeqIO
#from Bio.Seq import Seq

parser = argparse.ArgumentParser(prog='tbAnnotator.py', usage='%(prog)s [options]', description='description',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-gff', '--gff', type=str, default="annotation/H37RV_V5.gff3", help='Gff file with genes annotation')
parser.add_argument('vcf', type=str, help='Vcf file with snps')
parser.add_argument('-ddb', '--dreamDb', type=str, default="databases/tbDream.db", help='Tb dreamDb in hgvs-like format')
parser.add_argument('-tdb', '--tbvarDb', type=str, default="databases/tbVar.db", help='Tb tbVarDb')
parser.add_argument('-oj', '--outjson', type=str, default="result.json", help='Out json file with results')

args = parser.parse_args()
if not os.path.isfile(args.vcf):
	print('Bad vcf file')
	sys.exit(1)
if not os.path.isfile(args.dreamDb):
	print('Bad dreamDb file')
	sys.exit(1)
if not os.path.isfile(args.tbvarDb):
	print('Bad tbVarDb file')
	sys.exit(1)
if not os.path.isfile(args.gff):
	print('Bad gff file')
	sys.exit(1)


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
		if len(arr) and arr[2].strip() == 'gene':
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

def readDreamDb(fileName, genes):
	genes.sort(key=lambda x: x.start)
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


def readTbVarDb(fileName, genes):
	genes.sort(key=lambda x: x.start)
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
	#snvs.sort(key=lambda x: x.position)
	res = {}
	for snpPos in vcfSnps:
		filtered = []
		for snv in snvs:
			if snv.position == snpPos:
				filtered.append(snv)
		res[snpPos] = filtered
	return res

def getGeneDrugAssociation(genes, snvs):
	assoc = {}
	for snpId in snvs:
		geneId = snvs[snpId].gid
		drug = snvs[snpId].drug
		if geneId not in assoc:
			assoc[geneId] = {}
		if drug not in assoc[geneId]:
			assoc[geneId][drug] = 0
		assoc[geneId][drug] += 1
	for geneId in assoc:
		gene = getGene(genes, geneId)
		for drug in assoc[geneId]:
			assoc[geneId][drug] = assoc[geneId][drug] * 100. / abs(gene.start - gene.end)
	return assoc
def getGeneSeq(fasta, gene):
	print(len(fasta), gene.start, gene.end)
	subseq = fasta[gene.start: gene.end]
	if gene.strand < 0:
		subseq = subseq.reverse_complement()
	return subseq

def findGeneForSnp(genes, snp):
	for gene in genes:
		if gene.start <= snp.pos <= gene.end:
			return gene
	return None
def getDrugs(dreamSnvs):
	drugs = set()
	for snp in dreamSnvs:
		drugs.add(dreamSnvs[snp].drug)
	return list(drugs)



genes = processGffGenes(args.gff)
#referenceFasta = getReferenceFasta("H37RV_V5.fasta")
dreamSnvs = readDreamDb(args.dreamDb, genes)
tbvarSnvs = readTbVarDb(args.tbvarDb, genes)


assoc = getGeneDrugAssociation(genes, dreamSnvs)
allDrugs = getDrugs(dreamSnvs)
#print(allDrugs)

# for gene in assoc:
# 	print(gene)
# 	for drug in assoc[gene]:
# 		print('\t', drug, assoc[gene][drug])
def dumpJson(file, data):
	with open(file, 'w') as out:
		json.dump(data, out)
	return 0
def saveResults(file, results):
	dumpJson(file, results)


drugsScoreAssociated = {}
for drug in allDrugs:
	drugsScoreAssociated[drug] = 0

vcfSnps = readVcf(args.vcf)
coutLine = 1
resultJsonObj = {}
resultJsonObj['file'] = os.path.split(args.vcf)[1]
resultJsonObj['data'] = []
for snp in vcfSnps:
	drugs = []
	evidence = 0
	geneId = ''
	siftPred = ''
	freq = ''
	pubmed = ''

	gene = findGeneForSnp(genes, vcfSnps[snp])
	if gene:
		geneId = gene.gid
		if gene.gid in assoc:
			evidence = 1
			for drug in assoc[gene.gid]:
				drugsScoreAssociated[drug] = max(assoc[gene.gid][drug], drugsScoreAssociated[drug])
				drugs.append({'drug':drug, 'score': round(assoc[gene.gid][drug], 1)})
	if snp in tbvarSnvs:
		evidence = 2
		geneId = tbvarSnvs[snp].gene.gid
		freq = tbvarSnvs[snp].freq
		siftPred = tbvarSnvs[snp].siftPred
	if snp in dreamSnvs:
		evidence = 3
		geneId = dreamSnvs[snp].gid
		drugsScoreAssociated[drug] = max(100, drugsScoreAssociated[drug])
		drugs.append({'drug':dreamSnvs[snp].drug, 'score':100})
		pubmed =  dreamSnvs[snp].pubmed
	print(coutLine, evidence, vcfSnps[snp].pos, vcfSnps[snp].ref, vcfSnps[snp].alt, geneId, end = '\t', sep = '\t')
	if len(drugs):
		for drug in sorted(drugs, key=lambda x: -x['score']):
			print(drug['drug'], ':', round(drug['score'], 1), ';', sep = '', end = '')
	else:
		print('-', end ='')
	print('\t', siftPred, pubmed, sep = '\t')
	coutLine += 1
	resultJsonObj['data'].append({
		'position': vcfSnps[snp].pos,
		'refNucl': vcfSnps[snp].ref,
		'altNucl': vcfSnps[snp].alt,
		'geneId':geneId,
		'drugs':drugs,
		'evidence':evidence,
		'siftPred':siftPred,
		'freq':freq,
		'pubmed':pubmed
	})
resultJsonObj['data'].sort(key=lambda x: -x['evidence'])

resultJsonObj['drugsList'] = []
for drug in sorted(drugsScoreAssociated.keys(), key=lambda x: drugsScoreAssociated[x]):
	resultJsonObj['drugsList'].append({'drug':drug, 'resistanceScore':drugsScoreAssociated[drug]})


saveResults(args.outjson, resultJsonObj)
