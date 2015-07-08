import os,sys,getopt
from sys import path
path.append(sys.path[0]+"/lib/")
from FuncbioOTU import getAttriValueFromSeqName,fasta2dict

def fasta2DictDere(input_fasta,tag,seqdict):
    """ read a fasta into dictionary and remove the full equal(sequence and length)"""
    finput=open(input_fasta,"r")
    first=True
    for line in finput:
        if line.startswith(">"):
            if first:
                first=False
            else:
                name_list=seqdict.get(sequence,False)
                if name_list:
                    seqdict[sequence].append(name)
                else:
                    seqdict[sequence]=[name]
            #new circle
            name=line.strip().lstrip(">")+";SampleTag=%s;"%(tag)
            sequence=""
        else:
            line=line.strip()
            sequence+=line
    #the last
    name_list=seqdict.get(sequence,False)
    if name_list:
        seqdict[sequence].append(name)
    else:
        seqdict[sequence]=[name]
    return seqdict

def createLenDict(seqDict):
    """create a dictionary and sequence length as the key, a list that consisted by sequence and sequence name as value"""
    lenDict={}
    for sequence in seqDict:
        length=len(sequence)
        if lenDict.get(length,False):
            lenDict[length][sequence]=seqDict[sequence]
        else:
            lenDict[length]={sequence:seqDict[sequence]}
    return lenDict

def subseqRemove(dictA,dictB,Alen,Blen):
    
    #{len:{seq1:[name1,name2]}, seq2:[name3,name4]}  }
    "listA containing the longer sequences"
    #create A dict (key:seqname value:sequence)
    subDict={}
    for seqA in dictA:
        for n in range(Alen-Blen+2):
            subseq=seqA[n:(n+Blen)]
            subDict[subseq]=seqA
    #find subseq
    seqBlist=dictB.keys()
    for seqB in seqBlist:
        nameBList=dictB[seqB]
        if subDict.get(seqB,False):
            dictA[subDict[seqB]]+=nameBList
            del dictB[seqB]
        else:
            continue
    return dictA,dictB
def readFasta(filePathList):
    tagList=[]
    totalDict={}
    for line in open(filePathList,'r'):
        lineList=line.split()
        path=lineList[0]
        tag=lineList[1]
        tagList.append(tag)
        totalDict=fasta2DictDere(path,tag,totalDict)
    lenDict=createLenDict(totalDict)
    lenNum=len(lenDict)
    lenList=lenDict.keys()
    lenList.sort()
    lenList.reverse()
    a=0
    for n in range(lenNum):
        longList=lenDict[lenList[n]]
        for m in range(n+1,lenNum):
            shortList=lenDict[lenList[m]]
            dictAA,dictBB=subseqRemove(longList,shortList,lenList[n],lenList[m])
            lenDict[lenList[n]]=dictAA
            lenDict[lenList[m]]=dictBB
    return lenDict,tagList
if __name__=="__main__":
    usage="""usage: python dereplication.py 

    --list/-l   the list file (two column separated by tab) for all samples. The first column is the absolute path of sample files, and second column is the specified name of sample.
    --index/-n  to select the index ("samplesize" or "abundance") for dividing all unique sequences into both multiple and singleton tags.
    --help/-h help
    """
    character="sampleSize"
    opts,arg=getopt.getopt(sys.argv[1:],"l:n:h",['list=','index=','help'],)

    parameters=[a[0] for a in opts]
    if '-h' in parameters or '--help' in parameters:
        print usage
        sys.exit(1)
    if len(parameters)==0:
        print usage
        sys.exit(1)

    for i,a in opts:
        if i in ("--list","-l"):
            if not os.path.isfile(a):
                print "%s is not found."%(a)
                sys.exit(1)
            fastaList=a
        if i in ("--index","-n"):
            if a in ("samplesize","abundance"):
                character=a.strip()
                if character=="samplesize":
                    character="sampleSize"
            else:
                print "error, the value of %s must be 'sampleSize' or 'abundance'"%(i)
                sys.exit(1)
    outfasta="unique_sequences.fa"
    outtable="unique_sequences.table"
    outmultiple="unique_sequences_multiple.fa"
    outsingleton="unique_sequences_singleton.fa"
    unique_list="unique_sequences.list"

    fout=open(outfasta,'w')
    ftable=open(outtable,'w')
    fmulti=open(outmultiple,'w')
    fsingle=open(outsingleton,'w')
    flist=open(unique_list,'w')
    
    lenDict,tagList=readFasta(fastaList)
    ftable.write("sequenceName\tsampleSize\ttotalAbundance\taverageAbundance\t%s\n"%'\t'.join(tagList))
    n=0
    for lens in lenDict:
        for seq in lenDict[lens]:
            if len(lenDict[lens])==0:
                continue
            n+=1
            nameList=[]
            for seqname in lenDict[lens][seq]:
                sample=getAttriValueFromSeqName(seqname,'SampleTag')
                nameList.append(sample)
            distriList=[]
            for tag in tagList:
                count=nameList.count(tag)
                distriList.append(count)
            samplesize=len(tagList)-distriList.count(0)
            abundance=sum(distriList)
            average=float(sum(distriList))/len(distriList)
            flist.write(">unique%d;sampleSize=%d;abundance=%d;\t%s\n"%(n,samplesize,abundance,','.join(lenDict[lens][seq])))
            fout.write(">unique%d;sampleSize=%d;abundance=%d;\n%s\n"%(n,samplesize,abundance,seq))
            if character=="sampleSize":
                char=samplesize
            elif character=="abundance":
                char=abundance
            if char>1:
                fmulti.write(">unique%d;sampleSize=%d;abundance=%d;\n%s\n"%(n,samplesize,abundance,seq))
            else:
                fsingle.write(">unique%d;sampleSize=%d;abundance=%d;\n%s\n"%(n,samplesize,abundance,seq))
            ftable.write("unique%d\t%d\t%d\t%f\t%s\n"%(n,samplesize,abundance,average,'\t'.join([str(a) for a in distriList])))
    print "\noutput file:\n"
    print os.path.basename(outfasta)
    print os.path.basename(outtable)
    print os.path.basename(outmultiple)
    print os.path.basename(outsingleton)
    print os.path.basename(unique_list)
    print ''
