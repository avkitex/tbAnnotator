from __future__ import print_function
from Cheetah.Template import Template


def createHtmlReport(outFile, results, template = 'template.html'):
	templateDef=''.join(open(template).readlines())
	nameSpace = {'filename':results["file"], "data": results["data"], 'drugsList':results['drugsList']}
	t = Template(templateDef, searchList=[nameSpace])
	oHandle=open(outFile, 'w')
	print(t, file=oHandle)
	oHandle.close()
