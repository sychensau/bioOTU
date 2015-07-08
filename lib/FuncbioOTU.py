import sys,os
def fasta2dict(input_fasta):
    """convert a fasta type file into a python dictionary. sequence name as key, sequence as value."""
    if os.path.isfile(input_fasta):
        pass
    else:
        print "Error,%s not found."%(input_fasta)
        sys.exit(1)
    fasta_dict={}
    finput=open(input_fasta,"r")
    for line in finput:
        if line.startswith(">"):
            sequence_name=line.strip().lstrip(">")
            fasta_dict[sequence_name]=""
        else:
            sequence=line.strip()
            fasta_dict[sequence_name]+=sequence
    return fasta_dict
def getAttriValueFromSeqName(seqname,attribute):
    """get a specified attribute value from the sequence name."""
    attributeStart=seqname.find(attribute)
    if attributeStart==-1:
        return False
    else: 
        attributeEnd=attributeStart+len(attribute)
        valueStart=seqname.find("=",attributeEnd)+len("=")
        valueEnd=seqname.find(";",valueStart)
        if valueEnd>valueStart:
            value=seqname[valueStart:valueEnd].strip()
        else:
            value=seqname[valueStart:].strip()
        if len(value)==0:
            return False
        else:
            return value
def sortSequenceBySampleSize(inputfasta,outputfasta):
    if os.path.isfile(inputfasta):
        pass
    else:
        print "Error, %s not found."%(inputfasta)
        sys.exit(1)
    fout=open(outputfasta,"w")
    fastaDict=fasta2dict(input_fasta) #key:sampleSize, value:sampleNameList(which is the same sample size)
    sizeDict={}
    sequenceNameList=fastaDict.keys()
    for sequenceName in sequenceNameList:
        sampleSize=getAttriValueFromSeqName(sequenceName,"sampleSize")
        if sampleSize:
            pass
        else:
            print sequenceName,"\n"
            print "samplesSize not found."
            sys.exit(1)
        dictValue=sizeDict.get(sampleSize,None)
        if dictValue:
            sizeDict[sampleSize].append(sampleSize)
        else:
            sizeDict[sampleSize]=[sampleSize]
    sizeList=sizeDict.keys()
    sizeList.sort()
    sizeList.reverse()
    for size in sizeList:
        seqNameList=sizeDict[size]
        for seqName in seqNameList:
            sequence=fastaDict[seqName]
            fout.write(">%s\n"%(seqName))
            fout.write(sequence+"\n")
    return True
def sortOTUname(reStrList,NameList):
    newList=[]
    for reStr in reStrList:
        sortDict={}
        for line in NameList:
            if line.startswith(reStr):
                num=int(line.strip(reStr))
                sortDict[num]=line
        numList=sortDict.keys()
        numList.sort()
        for num in numList:
            newList.append(sortDict[num])
    return newList

def orderbySSAndAband(NameList):
    sizeVList=[]
    abundVList=[]
    valueDict={}
    value_nameDict={}
    reList=[]
    for Name in NameList:
        size=int(getAttriValueFromSeqName(Name,"sampleSize"))
        abundance=int(getAttriValueFromSeqName(Name,"abundance"))
        valueDict[Name]=[size,abundance]
        abundVList.append(abundance)
    abMax=max(abundVList)+1 #as weight value of sampleSize
    for name in valueDict:
        size=valueDict[name][0]
        abundance=valueDict[name][1]
        key=size*abMax+abundance
        if value_nameDict.get(key,False):
            value_nameDict[key].append(name)
        else:
            value_nameDict[key]=[name]

    valueList=value_nameDict.keys()
    valueList.sort()
    for value in valueList:
        reList+=value_nameDict[value]
    return reList
            
