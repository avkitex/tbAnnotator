import argparse

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

def readGffGenes(fileName):
	genes = []
	fHandle = open(fileName)
	rtype = 2
	gstart = 3
	gend = 4
	gstrand = 6
	gdescription = 8
	for line in fHandle:
		arr = line.strip().split()
		if len(arr):
			gid = ""
			for i in arr[gdescription].strip().split(';'):
				if "ID=" in i:
					gid = i.strip()[3:]
			if arr[gstrand].strip() == '+':
				strand = 1
			else:
				strand = -1
			genes.append(geneRecord(gid, int(arr[gstart].strip()), (arr[gend].strip()), strand, arr[gdescription]))
	fHandle.close()
	return genes

genes = readGffGenes(args.gff)
for i in genes:
	print(i.gid, i.start, i.end)
