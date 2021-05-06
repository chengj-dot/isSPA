#!/usr/bin/env python

import math, os, sys

try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(input1,num,meta,output) = parse_command_line()
        g=open(input1,"r")
	oo=open(output,"w")
        PI=3.14159265359
        inlst_line1=g.readlines()
        if((len(inlst_line1)-meta)<num):
             num=len(inlst_line1)-meta
        
        score=[([0]*2) for i in range(len(inlst_line1)-meta)]
        lst_line=len(inlst_line1[4].split())
        for j in range(2,(lst_line)):
            tmp=str(inlst_line1[4].split('\t')[j].split('=')[0])
            if(tmp=="score"):
	        index=int(j)
        for i in range(meta,len(inlst_line1)):
            nn=int(i)
            score[i-meta][0]=nn
            score[i-meta][1]=abs(float(inlst_line1[i].split('\t')[index].split('=')[1]))
        score_sorted=sorted(score,key=lambda x:x[1],reverse=True);
        for j in range(0,num):
            nn=score_sorted[j][0]
            oo.write(inlst_line1[nn])
        g.close()
        oo.close()
		
def parse_command_line():
	usage="%prog <input1> <number> <metadata line> <output>"
	parser = OptionParser(usage=usage, version="%2016, by Kyouko>")
	
	if len(sys.argv)<3: 
		print "<input1> <number> <metadata line> <output> "
		sys.exit(-1)
	
	(options, args)=parser.parse_args()
        input1=args[0]
        num=int(args[1])
        meta=int(args[2])
	output=args[3]
	return (input1,num,meta,output)


if __name__== "__main__":
	main()


