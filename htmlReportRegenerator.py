import sys, argparse, os, json
from htmlGenerator import createHtmlReport

parser = argparse.ArgumentParser(prog='htmlReportRegenerator.py', usage='%(prog)s [options]', description='Regeneration of html report',
								 epilog="\xa9 Avktex 2016")
parser.add_argument('-inj', '--inputJson', type=str, help='Input saved fson info', required=True)
parser.add_argument('-htmlrep', '--htmlreport', type=str, default='report.html',  help='Html out file.')
parser.add_argument('-tp', '--template', type=str, default='template.html', help='Path to template')

args = parser.parse_args()
def loadJson(file):
	try:
		with open(file, 'r') as inp:
			data = json.load(inp)
		return data
	except Exception as e:
		print("Not able to read json file. " + str(e))
		return None

if not os.path.isfile(args.template):
	print('Not able to locate template file')
	sys.exit(1)
if not os.path.isfile(args.inputJson):
	print('Not able to locate input json file')
	sys.exit(1)

results=loadJson(args.inputJson)

createHtmlReport(args.htmlreport, results, args.template)
