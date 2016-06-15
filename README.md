# tbAnnotator

A programm complex for annotating mycobacterium tuberculosis SNPs with information about drug resistance. 
It also performes scoring of SNPs not found in provided databases

Developed by avkitex, phill.gusev, pankevich-ev on Genehack-2


## Requirenments
```sh
apt-get install python2.7, python-pip
pip install Cheetah
```
## RUN

```sh
python tbAnnotator.py patient.vcf [ -gff annotation_gff -ddb dream_TB_database.db -tdb tbvar_database.db -oj output.json ]
```
output.json will contain final output information whick can be visualized with in html.
```sh
python htmlReportRegenerator.py -inj output.json -htmlrep report.html
```
report.html will contain sorted table with all drugs colored according to predicted drug resistance and a table with all snps and all information about them.
