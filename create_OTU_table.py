import sys
import os
import getopt
import random
import commands
import multiprocessing
from sys import path
path.append(sys.path[0]+"/lib/")
from FuncbioOTU import getAttriValueFromSeqName,fasta2dict,sortOTUname

def table2dict(tablefile):
    reDict={}
    for line in open(tablefile,"r"):
        lineList=line.strip().split()
        key=lineList[0]
        value=lineList[1:]
        reDict[key]=value
    return reDict

def sumbyotu(tableDict,otuFile):
    newDict={}
    for line in open(otuFile,'r'):
        lineList=line.strip().strip(',').split()
        OTUname=lineList[0]
        OTUseqList=lineList[1].split(',')
        newDict["sequenceName"]=tableDict["sequenceName"]
        sumList=[]
        first=True
        for seqName in OTUseqList:
            i=seqName.index(';')
            seqName=seqName[:i]
            numList=tableDict[seqName]
            if first:
                sumList=[float(a) for a in numList]
                first=False
            else:
                for n in range(len(sumList)):
                    sumList[n]+=float(numList[n])
        newDict[OTUname]=sumList

    #count sample size
    for otu in newDict:
        if otu=="sequenceName":
           pass
        else:
            groupSize=0
            nonSize=0
            k=0
            for i in newDict[otu][3:]:
                if i!=0:
                    groupSize+=1
                else:
                    nonSize+=1
                newDict[otu][3+k]=int(i)
                k+=1
            average=float(newDict[otu][1])/float(groupSize+nonSize)
            newDict[otu][0]=groupSize
            newDict[otu][2]=average
            newDict[otu][1]=int(newDict[otu][1])
    return newDict

def sumbygroup(tableDict,groupFile):
    groupDict={}
    newDict={}
    for line in open(groupFile,'r'):
        lineList=line.strip().split()
        if groupDict.get(lineList[1],False):
            groupDict[lineList[1]].append(lineList[0])
        else:
            groupDict[lineList[1]]=[lineList[0]]
    First=True
    for key in tableDict:
        valueList=tableDict[key]
        newDict[key]=valueList[0:3]
        for group in groupDict:
            groupvalueList=[]
            if key=="sequenceName":
                newDict["sequenceName"].append(group)
            else:
                for sampleName in groupDict[group]:
                    indexId=tableDict["sequenceName"].index(sampleName)
                    groupvalueList.append(float(valueList[indexId]))
                total=sum(groupvalueList)
                newDict[key].append(total)
    #count group size
    for otu in newDict:
        if otu=="sequenceName":
           newDict[otu][0]="groupSize"
        else:
            groupSize=0
            nonSize=0
            k=0
            for i in newDict[otu][3:]:
                if i!=0:
                    groupSize+=1
                else:
                    nonSize+=1
                newDict[otu][3+k]=int(i)
                k+=1
            average=float(newDict[otu][1])/float(groupSize+nonSize)
            newDict[otu][0]=groupSize
            newDict[otu][2]=average
            newDict[otu][1]=int(newDict[otu][1])
    return newDict

if "__main__"==__name__:
    usage='''
    --OTUlist/-u   list files (file names are separated by ',') of these constructed OTUs by both taxonomy-guided and heuristic search methods ("taxonomy_guided_OTU.list,heuristic_search_OTU.list")

    --table/-t   table file containing information of both sample size and abundance for all unique sequences created by module of "dereplication.py" ("unique_sequences.table")

    --group/-g   list file provided by user, which contains tab-separated two column: the first column is the specified name of sample, and the second column is group ID

'''
    knowOTU=False
    unknowOTU=False
    group=False
    opts,arg=getopt.getopt(sys.argv[1:],"u:t:g:h",["OTUlist=","table=","group=","help"])

    parameters=[a[0] for a in opts]
    if '-h' in parameters or '--help' in parameters:
        print usage
        sys.exit(1)
    if len(parameters)==0:
        print usage
        sys.exit(1)
    if '-u' not in parameters and '--OTUlist' not in parameters:
        print "***Error, an OTU list file is requred.***\n"
        print usage
        sys.exit(1)
    if '-t' not in parameters and '--table' not in parameters:
        print "***Error, a table file is requred.***\n"
        print usage
        sys.exit(1)

    for i,a in opts:
        if i in ('-u','--OTUlist'):
            OTUfileList=a.strip().split(',')
            for OTUfile in OTUfileList:
                if not os.path.isfile(OTUfile):
                    print "%s is not found."%(OTUfile)
                    sys.exit(1)
        if i in ('-t','--table'):
            if not os.path.isfile(a):
                print "%s is not found."%(a)
                sys.exit(1)
            else:
                table=a
        if i in ('-g','--group'):
            if not os.path.isfile(a):
                print "%s is not found."%(a)
                sys.exit(1)
            else:
                group=a

    sample_output="sample_OTU.table"
    fout=open(sample_output,"w")

    tableDict=table2dict(table)

    #by sample
    table1={}
    for OTUfile in OTUfileList:
        resultTable=sumbyotu(tableDict,OTUfile)
        table1.update(resultTable)
    otuNameList=table1.keys()
    otuNameList.sort()

    #out result
    otuNameList.remove("sequenceName")
    head=table1["sequenceName"]

    fout.write("%s\t%s\n"%("sequenceName","\t".join(head)))
    otuNameList=sortOTUname(['OTU','hsOTU'],otuNameList)
    for otu in otuNameList:
        otuCount=table1[otu]
        fout.write(otu+"\t")
        for num in otuCount:
            fout.write(str(num)+"\t")
        fout.write("\n")

    #by group
    if group:
        table1=sumbygroup(table1,group)
        group_output="group_OTU.table"
        fout=open(group_output,"w")
        otuNameList=table1.keys()
        otuNameList.sort()

    #out result
        otuNameList.remove("sequenceName")
        head=table1["sequenceName"]

        fout.write("%s\t%s\n"%("sequenceName","\t".join(head)))
        otuNameList=sortOTUname(['OTU','hsOTU'],otuNameList)
        for otu in otuNameList:
            print otu
            otuCount=table1[otu]
            fout.write(otu+"\t")
            for num in otuCount:
                fout.write(str(num)+"\t")
            fout.write("\n")
    print "\noutput file:\n"
    print os.path.basename(sample_output)
    if group:
        print os.path.basename(group_output)
    print ''
