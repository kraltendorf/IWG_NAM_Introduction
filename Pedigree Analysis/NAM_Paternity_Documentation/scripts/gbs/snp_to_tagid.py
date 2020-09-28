#!/usr/bin/env python
import argparse
import re

#From Liang Gao 

parser = argparse.ArgumentParser (description = "SNP id to tag sequence retrieval after TASSEL GBSv2")
parser.add_argument ('-t','--tagseq', help ="tag sequences extracted using: perl ~/tassel/tassel-5-standalone/run_pipeline.pl -fork1 -GetTagSequenceFromDBPlugin -db YOURSQLITE.db -o tagseq.txt -endPlugin -runfork1" )
parser.add_argument ('-s', '--sql_table_export', help = 'SQL table join using: sqlite3 -separator $\'\\t\' YOURSQLITE.db < /homes/lianggao/scripts/SQL/GBSv2_table_join_to_get_tagid.sql > YOURPrefix_sql_join.txt')
args = parser.parse_args()

# tagseq = "/homes/lianggao/gbs/projects/LakinFuller/results_ref1.0_end_to_end_very_sens/tagseq.txt"
# sql_table_export = "/homes/lianggao/gbs/projects/LakinFuller/results_ref1.0_end_to_end_very_sens/LakinFuller_sql_join.txt"

tagseq = args.tagseq
sql_table_export = args.sql_table_export


#########(1) load tagid to seq file get a hash
### seq[tagid]="atcg..."

fh = open (tagseq,'rU')
header = fh.readline()
tagid = 1
seq = dict()
for line in fh:
    line = line.rstrip()
    seq[tagid]=line
    tagid=tagid+1
fh.close()
    

#(2) read in the joined sql table

fh = open (sql_table_export,'rU')

for line in fh:
    line = line.rstrip()
    F = line.split()
    (chr,pos,tagid)=(F[1],F[2],F[5])
    ID = chr+"_"+pos

    newID =">chr"+ID+"_tag"+tagid
    print (newID)
    print (seq[int(tagid)])

