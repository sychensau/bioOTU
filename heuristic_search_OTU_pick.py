import sys
import os
import getopt
import stat
import time
import random
import commands
import multiprocessing
from sys import path
path.append(sys.path[0]+"/lib/")
from FuncbioOTU import getAttriValueFromSeqName,fasta2dict,orderbySSAndAband

def multi(exe_path,fastafile,output_file_path,filerDistance):
    mycommand=exe_path+" %s %s %f" %(fastafile,output_file_path,filerDistance)
    commands.getoutput(mycommand)

def dist2Dict(distFile):
    reDict={}
    for line in open(distFile,"r"):
        lineList=line.strip().split()
        dist=lineList[2]
        if float(dist)==0:
            pass
        else:
            name1=lineList[0]+lineList[1]
            name2=lineList[1]+lineList[0]
            reDict[name1]=float(dist)
            reDict[name2]=float(dist)
    return reDict

def listbySampleSize(fastafile):
    reDict={}
    for line in open(fastafile,"r"):
        if line.startswith(">"):
            size=getAttriValueFromSeqName(line,"sampleSize")
            if reDict.get(size,False):
                reDict[size].append(line.lstrip(">").rstrip())
            else:
                reDict[size]=[line.lstrip(">").rstrip()]
    return reDict

def cluster(nameListSorted,distDict,threshold):
    reDict={}
    lenList=len(nameListSorted)
    for n in range(lenList):
        minDist=1.1
        minSeqName=""
        for m in range(n+1,lenList):
            nameKey=nameListSorted[n]+nameListSorted[m]
            distValue=distDict.get(nameKey,1)
            if distValue<=minDist:
                minDist=distValue
                minSeqName=nameListSorted[m]
        if minDist<=threshold:
            subSeqList=reDict.get(minSeqName,[])
            subSeqList.append(nameListSorted[n])
            if reDict.get(nameListSorted[n],False):
                subSeqList+=reDict[nameListSorted[n]]
                del reDict[nameListSorted[n]]
            else:
                pass
            reDict[minSeqName]=subSeqList
        else:
            if reDict.get(nameListSorted[n],False):
                pass
            else:
                reDict[nameListSorted[n]]=[]
    return reDict

def sequenceFilter(seqDict,minSize,minAbundance):
    nameList=seqDict.keys()
    for name in nameList:
        size=int(getAttriValueFromSeqName(name,"sampleSize"))
        abundance=int(getAttriValueFromSeqName(name,"abundance"))
        if size<minSize or abundance<minAbundance:
            del seqDict[name]
    return seqDict

def multiForOneGenus_1(exe_path,seqfile,output_distance,filerDistance):
    mycommand=exe_path+" %s %s %f" %(seqfile,output_distance,filerDistance)
    commands.getoutput(mycommand)

def multiForOneGenus_2(exe_path,afile,bfile,output_file_path,filerDistance):
    usearch_command=exe_path+" %s %s %s %f" %(afile,bfile,output_file_path,filerDistance)
    commands.getoutput(usearch_command)

def multprocessdistance(fastaDict,tempDir,processors,script_loc,filerDistance):

    #sequence filtering by sample size and abundance
    #calculate distance start
    dictLen=len(fastaDict)
    keyList=fastaDict.keys()
    oneLen=dictLen/processors+1
    k=1
    one=0
    subfileList=[]
    while len(keyList)>0:
        subfile=tempDir+"/inputseed_"+str(k)+".fasta"
        subfileList.append(subfile)
        fout=open(subfile,"w")
        k+=1
        for n in range(oneLen):
            if len(keyList)>0:
                name=keyList.pop()
                seq=fastaDict[name]
                fout.write(">%s\n%s\n"%(name,seq))
            else:
                break
        fout.close()

    #distance calculate in one file
    resultList=[] #store distance file path
    distanceEXE=script_loc+"/lib/distanceCalculateWithKmerOneFile"
    pool=multiprocessing.Pool(processes=processors)
    for fastafilepath in subfileList:
        distanceFile=tempDir+"/%s.distance"%(os.path.basename(fastafilepath))
        resultList.append(distanceFile)
        pool.apply_async(multiForOneGenus_1,(distanceEXE,fastafilepath,distanceFile,filerDistance))
    pool.close()
    pool.join()
    #distance calculate bewteen two file
    distanceEXE=script_loc+"/lib/distanceCalculateWithKmer"
    pool=multiprocessing.Pool(processes=processors)
    for b in range(len(subfileList)):
        for c in range(b+1,len(subfileList)):
            distanceFile=tempDir+"/%d_%d.distance"%(b,c)
            resultList.append(distanceFile)
            before_file=subfileList[b]
            next_file=subfileList[c]
            pool.apply_async(multiForOneGenus_2,(distanceEXE,before_file,next_file,distanceFile,filerDistance))
    pool.close()
    pool.join()
    #combine result 
    distanceResult=tempDir+"/unknown.distance"
    fout=open(distanceResult,"w")
    for resultfile in resultList:
        for line in open(resultfile,'r'):
            fout.write(line)
    fout.close()
    return distanceResult

def multi(exe_path,seedfile,candidatefile,output_file_path,filerDistance):
    mycommand=exe_path+" %s %s %s %f" %(seedfile,candidatefile,output_file_path,filerDistance)
    commands.getoutput(mycommand)

def distanceCalculate(seed,candidate,filerDistance,proc,tempDir,script_loc):

    output=tempDir+"/seed_candidate.distance"

    exepath=script_loc+"/lib/distanceCalculateWithKmer"

    fastaDict=fasta2dict(seed)
    dictLen=len(fastaDict)
    keyList=fastaDict.keys()
    oneLen=dictLen/proc+1
    k=1
    one=0
    subSeedPathList=[]
    while len(keyList)>0:
        subfile=tempDir+"/inputseed_"+str(k)+".fasta"
        subSeedPathList.append(subfile)
        fout=open(subfile,"w")
        k+=1
        for n in range(oneLen):
            if len(keyList)>0:
                name=keyList.pop()
                seq=fastaDict[name]
                fout.write(">%s\n%s\n"%(name,seq))
            else:
                break
        fout.close()
    outputfileList=[]
    pool=multiprocessing.Pool(processes=proc)
    k=0
    for subfile in subSeedPathList:
        k+=1
        outputfile=tempDir+"/outdistance_"+str(k)+".distance"
        outputfileList.append(outputfile)
        pool.apply_async(multi,(exepath,subfile,candidate,outputfile,filerDistance))
    pool.close()
    pool.join()
    fout=open(output,'w')
    for suboutput in outputfileList:
        for line in open(suboutput,"r"):
            fout.write(line)
    fout.close()
    return output

def assignedSingleton2multiple(seedFile,distanceFile,threshold,singletonfile):
    singletonList=[] #return 
    reDict={}
    usesingleton=[]
    seedDict={}
    for line in open(seedFile,"r"):
        if line.startswith(">"):
            seedName=line.strip().strip(">")
        # 
        if seedDict.get(seedName,"a")==[]:
            pass
        else:
            seedDict[seedName]=[]

    for line in open(distanceFile,'r'):
        lineList=line.strip().split()
        candidateDict=reDict.get(lineList[1],None)
        if candidateDict:
            if candidateDict.get(lineList[2],None):
                reDict[lineList[1]][lineList[2]].append(lineList[0])
            else:
                reDict[lineList[1]][lineList[2]]=[lineList[0]]
        else:
            reDict[lineList[1]]={}
            reDict[lineList[1]][lineList[2]]=[lineList[0]]
    for Candi in reDict:
        candidateName=Candi
        candiDict=reDict[Candi]
        distanceList=candiDict.keys()
        minDistance=min(distanceList)
        if float(minDistance)<=float(threshold):
            usesingleton.append(Candi)
            seedList=candiDict[minDistance]
            if len(seedList)>1:
                maxSize=0
                selectSeedSequenceName=""
                for seqName in seedList:
                    sampleSize=getAttriValueFromSeqName(seqName,"sampleSize")
                    if int(sampleSize)>int(maxSize):
                        maxSize=sampleSize
                        selectSeedSequenceName=seqName
                    else:
                        pass
            else:
                selectSeedSequenceName=seedList[0]
            reDict[Candi]=selectSeedSequenceName
        else:
            reDict[Candi]='None'
    for line in open(singletonfile,'r'):
        if line.startswith('>'):
            name=line[1:].strip()
            if name in usesingleton:
                continue
            else:
                singletonList.append(name)
    #
    for candiName in reDict:
        ParentsSeq=reDict[candiName]
        if ParentsSeq=='None':
            pass
        else:
            seedDict[ParentsSeq].append(candiName)
    return seedDict,singletonList
if __name__=="__main__":
    usage='''
    --pengding_nonchimera/-s   (required) input file containing all pengding multiple sequences after chimera detection ("pengding_sequences_multiple.nonchimera")

    --pending_singleton/-g   (required) input file containing all pengding singleton sequences ("pending_sequence_singleton.fa")

    --threshold/-r   (optional) to specify the distance cutoff for OTUs clustering, default:0.03

    --kdistance/-k   (optional) to specify cutoff of kmer distance, default:0.5

    --processors/-p   (optional) processors, default:4
'''
    threshold=0.03
    filerDistance=0.5
    proc=4
    opts,arg=getopt.getopt(sys.argv[1:],"s:k:p:g:r:h",["pengding_nonchimera=","kdistance=","processors=","pending_singleton=","threshold=",'help'])

    parameters=[a[0] for a in opts]
    if '-h' in parameters or '--help' in parameters:
        print usage
        sys.exit(1)
    if len(parameters)==0:
        print usage
        sys.exit(1)
    if '-s' not in parameters and '--pengding_nonchimera' not in parameters:
        print "***Error, --assigned/-a is requred.***\n"
        print usage
        sys.exit(1)
    if '-g' not in parameters and '--pending_singleton' not in parameters:
        print "***Error, --unassigned/-u is requred.***\n"
        print usage
        sys.exit(1)

    for i,a in opts:
        if i in ('-s','--pengding_nonchimera'):
            if not os.path.isfile(a):
                print "***%s is not found.***"%(a)
                sys.exit(1)
            seqfile=a
        if i in ('-g','--pending_singleton'):
            if not os.path.isfile(a):
                print "***%s is not found.***"%(a)
                sys.exit(1)
            singleton=a
        if i in ("--threshold","-r"):
            try:
                threshold=float(a)
            except:
                print "***Error, the value of threshold (--threshold/-r) must be decimals.***\n"
                print usage
                sys.exit(1)
        if i in ("--kdistance","-k"):
            try:
                filerDistance=float(a)
            except:
                print "***Error, the kmer distance (--kdistance/-k) must be decimals.***\n"
                print usage
                sys.exit(1)
        if i in ('-p','--processors'):
            try:
                processors=int(a)
            except:
                print "***Error, the processors (--processors/-p) must be integer.***\n"
                print usage
                sys.exit(1)

    otulist="heuristic_search_OTU.list"
    outseq="heuristic_search_OTU.fa"
    script_loc=os.path.split(os.path.realpath(sys.argv[0]))[0]
    exepath=script_loc+"/lib/distanceCalculateWithKmerOneFile"
    tempDir=script_loc+"/temp"+str(random.randint(10000,99999))
    os.mkdir(tempDir)
    fastaDict=fasta2dict(seqfile)
    singleDict=fasta2dict(singleton)

    #Executable permissions
    muthurPath=script_loc+"/lib/Mothur.cen_64/mothur/mothur"
    distance1=script_loc+"/lib/distanceCalculateWithKmer"
    distance2=script_loc+"/lib/distanceCalculateWithKmerOneFile"

    try:
        os.chmod(muthurPath,stat.S_IRWXU)
    except:
        commands="chmod a+x %s"%(muthurPath)
        print "Please give executable permission to %s by this terminal commands:\n%s"%(muthurPath,commands)
        sys.exit(1)
    try:
        os.chmod(distance1,stat.S_IRWXU)
    except:
        commands="chmod a+x %s"%(distance1)
        print "Please give executable permission to %s by this terminal commands:\n%s"%(distance1,commands)
        sys.exit(1)
    try:
        os.chmod(distance2,stat.S_IRWXU)
    except:
        commands="chmod a+x %s"%(distance2)
        print "Please give executable permission to %s by this terminal commands:\n%s"%(distance2,commands)
        sys.exit(1)


    #calculate distance of multiple sequences and singleton sequences
    multipleSingletonDistance=distanceCalculate(seqfile,singleton,filerDistance,proc,tempDir,script_loc)
    seedDict,discardSingletonList=assignedSingleton2multiple(seqfile,multipleSingletonDistance,threshold,singleton)
    print "discard sequences(singleton):",len(discardSingletonList)
    #distance of multiple in onefile
    distanceAll=multprocessdistance(fastaDict,tempDir,proc,script_loc,filerDistance)
    distanceDict=dist2Dict(distanceAll)
    nameList=fastaDict.keys()
    nameListSorted=orderbySSAndAband(nameList)
    clusterDict=cluster(nameListSorted,distanceDict,threshold)

    flist=open(otulist,"w")
    k=0
    RE="hsOTU"
    fout=open(outseq,"w")
    for rep in clusterDict:
        k+=1
        multi_rep_list=clusterDict[rep]
        OTUlistStrList=[]
        OTUlistStrList.append("%s%d\t%s*"%(RE,k,rep))
        fout.write(">%s\n%s\n"%(rep,fastaDict[rep]))
        if seedDict.get(rep,False):
            if len(seedDict[rep])>0:
                OTUlistStrList.append(','.join(seedDict[rep]))
                for singleRepName in seedDict[rep]:
                    fout.write(">%s\n%s\n"%(singleRepName,singleDict[singleRepName]))
        for submulti in multi_rep_list:
            OTUlistStrList.append('%s'%(submulti))
            fout.write(">%s\n%s\n"%(submulti,fastaDict[submulti]))
            if seedDict.get(submulti,False):
                OTUlistStrList.append(','.join(seedDict[submulti]))
                for singlesub in seedDict[submulti]:
                    fout.write(">%s\n%s\n"%(singlesub,singleDict[singlesub]))
        flist.write(",".join(OTUlistStrList))
        flist.write('\n')
    flist.close()
    fout.close()
    remove_intermediate_file=r"rm -rf "+script_loc+r"/temp* "
    commands.getoutput(remove_intermediate_file)
    print "\noutput file:\n"
    print os.path.basename(otulist)
    print os.path.basename(outseq)
    print ''
