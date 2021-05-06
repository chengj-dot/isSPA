#!/usr/bin/env python

import math, os, sys, random
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(inlst1, inlst2, scale, output) =  parse_command_line()
	g = open(inlst1, "r")
        f = open(inlst2, "r")
	num=0
	oo=open(output,"w")
	inlst_line1=g.readlines()
        inlst_line2=f.readlines()
        k=3
	for i in range(0,len(inlst_line1)):
               # num = int(inlst_line1[i].split('\t')[0])
                m = int(inlst_line1[i].split('\t')[1].split('.')[1])
                #m=9
                if(1):
                    nn = inlst_line2[m+3].split('\t')[0]
                    s1 = inlst_line2[m+3].split('\t')[1]
                    s2 = inlst_line2[m+3].split('\t')[2]
                    s3 = inlst_line2[m+3].split('\t')[3]
                    s4 = float(inlst_line2[m+3].split('\t')[4].split('=')[1])
                    euler1 = float(inlst_line1[i].split('\t')[5].split('=')[1].split(',')[0])
                    euler2 = float(inlst_line1[i].split('\t')[5].split('=')[1].split(',')[1])
                    euler3 = float(inlst_line1[i].split('\t')[5].split('=')[1].split(',')[2])
                    score = float(inlst_line1[i].split('\t')[7].split('=')[1])
             #   m[i-3] = int(inlst_line1[i].split('\t')[0])
                    s5 = "euler="
                    s6 = "center="
                    cx = (float(inlst_line1[i].split('\t')[6].split('=')[1].split(',')[0]))*scale
                    cy = (float(inlst_line1[i].split('\t')[6].split('=')[1].split(',')[1]))*scale
                    #cx = float(inlst_line1[i].split('\t')[5])
                    #cy = float(inlst_line1[i].split('\t')[6])
                    oo.write(str(nn)+"\t"+str(s1)+"\t"+str(s2)+"\t"+str(s3)+"\t"+"dfang="+str(s4)+"\t"+str(s5)+str(euler1)+","+str(euler2)+","+str(euler3)+"\t"+str(s6)+str(cx)+","+str(cy)+"\t"+"score="+str(score)+"\n")
                    num += 1                           
        g.close()
        f.close()
	oo.close()

	
def parse_command_line():
	usage="%prog <input lstfile> <original particle file> <scale factor> <output>"
	parser = OptionParser(usage=usage, version="%1")
	
	if len(sys.argv)<4: 
		print "<input lstfile> <original particle file> <scale factor> <output>"
		sys.exit(-1)
	
	(options, args)=parser.parse_args()
	
	inlst1 = args[0]
	inlst2 = args[1]
        scale = int(args[2])
        output = args[3]
	return (inlst1, inlst2, scale, output)

if __name__== "__main__":
	main()
