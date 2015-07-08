import sys
import os
import getopt
import random
import commands
import stat
import multiprocessing
from sys import path
path.append(sys.path[0]+"/lib/")
from FuncbioOTU import getAttriValueFromSeqName,fasta2dict

##########################################################
# 1.calculate distance of annotated and unannotated.     #
# 2.assign unannotated sequence to annotated.            #
# 3.assign singleton sequences to annotated sequences.   #
# 4.pick OTU with annotated sequence. (singleton sequence#
#    and unannotated sequence which was assigned only as #
#     a repeat sequence of annotated.)                   #
#********************************************************#
# created in 2015.6.27                .                  #
#                                                        #
###########################################################

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
            reDict[name1]=dist
            reDict[name2]=dist
    return reDict

def clustSeedCandidate(seedfile,distanceFile,threshold,candidatefile):
    seedDict={}
    candidateList=[] #return 
    reDict={}
    useCandidate=[]
    for line in open(seedfile,"r"):
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
            useCandidate.append(Candi)
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
    for line in open(candidatefile,'r'):
        if line.startswith('>'):
            name=line[1:].strip()
            if name in useCandidate:
                continue
            else:
                candidateList.append(name)
    #
    for candiName in reDict:
        ParentsSeq=reDict[candiName]
        if ParentsSeq=='None':
            pass
        else:
            seedDict[ParentsSeq].append(candiName)
    return seedDict,candidateList

def createSeedTaxonomyDict(taxonomyFile,seedDict):
    returnDict={}
    taoxnomyDict={}
    for line in open(taxonomyFile,"r"):
        lineList=line.strip().split()
        taoxnomyDict[lineList[0]]=lineList[1]
    for seedName in seedDict:
        seedTaxonomy=taoxnomyDict[seedName]
        if returnDict.get(seedTaxonomy,None):
            returnDict[seedTaxonomy].append(seedName)
        else:
            returnDict[seedTaxonomy]=[seedName]
    return returnDict

def generateCountFile(inputSeqFilePath,countFilePaht):
    """used by cluster(mothur)"""
    fout=open(countFilePaht,'w')
    fout.write("%s\t%s\n"%("seqName","count"))
    for line in open(inputSeqFilePath,'r'):
        if line.startswith('>'):
            name=line.rstrip().lstrip('>')
            fout.write("%s\t%d\n"%(name,1))
    fout.close()
    return True

def generateCountFile2(seqNameList,countFilePath):
    """used by cluster(mothur)"""
    fout=open(countFilePath,'w')
    fout.write("%s\t%s\n"%("seqName","count"))
    for line in seqNameList:
        name=line.strip()
        fout.write("%s\t%d\n"%(name,1))
    fout.close()
    return True

def OTUcallIngenus(inputSeqFilePath,Distexepath,muthurPath):
    distanceFile=inputSeqFilePath+".distance"
    #create a count file for 
    countFilePath=inputSeqFilePath+".count"
    generateCountFile(inputSeqFilePath,countFilePath)
    #
    mycommand=Distexepath+" %s %s %f" %(inputSeqFilePath,distanceFile,1.1)
    commands.getoutput(mycommand)
    #cluster
    muthurParaStr="\"#cluster(column=%s, method=average, count=%s)\""%(distanceFile,countFilePath)
    returnMsg=commands.getoutput(muthurPath+" "+muthurParaStr)
    #deal with result
    return True

def outputRetainedCandidate(sourceFasta,candidateList,candidateSeq):
    fcand=open(candidateSeq,"w")
    sourceDict=fasta2dict(sourceFasta)
    for seqName in candidateList:
        sequence=sourceDict.get(seqName,False)
        if sequence:
            fcand.write(">"+seqName+"\n"+sequence+"\n")
        else:
            print "warnning, %s is not found in %s."%(seqName,sourceFasta)#there is some error.
    fcand.close()
    return True

def mothurOTUlist2dict(otulist):
    OTUlistdict={}
    for line in open(otulist,'r'):
        lineList=line.strip().split()
        key=lineList[0]
        if key=="unique":
            key=0
        value=lineList[2:]
        OTUlistdict[key]=value
    return OTUlistdict

def multiOTUcall(fastaFile,genusList,distanceEXE,muthurPath):
    fout=open(fastaFile,"w")
    for seqName in genusList:
        sequence=seedSeqDict[seqName]
        fout.write(">%s\n%s\n"%(seqName,sequence))
    fout.close()
    OTUcallIngenus(fastaFile,distanceEXE,muthurPath)

def multi(exe_path,seedfile,candidatefile,output_file_path,filerDistance):
    usearch_command=exe_path+" %s %s %s %f" %(seedfile,candidatefile,output_file_path,filerDistance)
    commands.getoutput(usearch_command)

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

def multiForOneGenus_1(exe_path,seqfile,output_distance):
    filerDistance=1.1 #not filter
    mycommand=exe_path+" %s %s %f" %(seqfile,output_distance,filerDistance)
    commands.getoutput(mycommand)

def multiForOneGenus_2(exe_path,afile,bfile,output_file_path):
    filerDistance=1.1 #not filter
    usearch_command=exe_path+" %s %s %s %f" %(afile,bfile,output_file_path,filerDistance)
    commands.getoutput(usearch_command)

def multprocessForOneGenus(seqDict,oneGenus,distanceEXE,muthurPath,processors,k,tempDir,script_loc):
    seqNameList=list(oneGenus)
    seqNum=len(seqNameList)
    seqNumPerProcess=int(seqNum/processors)+1
    n=0
    subfileList=[]
    while len(seqNameList)>0:
        n+=1
        subfile=tempDir+"/subdistance_%d_%d.fa"%(k,n)
        subfileList.append(subfile)
        fout=open(subfile,"w")
        for m in range(seqNumPerProcess):
            if len(seqNameList)>0:
                seqName=seqNameList.pop()
                seq=seqDict[seqName]
                fout.write(">%s\n%s\n"%(seqName,seq))
            else:
                break
        fout.close()
    #distance calculate in one file
    resultList=[] #store distance file path
    distanceEXE=script_loc+"/lib/distanceCalculateWithKmerOneFile"
    a=0
    pool=multiprocessing.Pool(processes=processors)
    for fastafilepath in subfileList:
        a+=1
        distanceFile=tempDir+"/OneGenus_%d_%d.fa.distance"%(k,a)
        resultList.append(distanceFile)
        pool.apply_async(multiForOneGenus_1,(distanceEXE,fastafilepath,distanceFile))
    pool.close()
    pool.join()
    #distance calculate bewteen two file
    distanceEXE=script_loc+"/lib/distanceCalculateWithKmer"
    pool=multiprocessing.Pool(processes=processors)
    for b in range(len(subfileList)):
        for c in range(b+1,len(subfileList)):
            a+=1
            distanceFile=tempDir+"/OneGenus_%d_%d.fa.distance"%(k,a)
            resultList.append(distanceFile)
            before_file=subfileList[b]
            next_file=subfileList[c]
            pool.apply_async(multiForOneGenus_2,(distanceEXE,before_file,next_file,distanceFile))
    pool.close()
    pool.join()
    #combine result 
    distanceResult=tempDir+"/OneGenus_%d.fa.distance"%(k)
    fout=open(distanceResult,"w")
    for resultfile in resultList:
        for line in open(resultfile,'r'):
            fout.write(line)
    fout.close()
    #create a count file for 
    countFilePath=distanceResult+".count"
    generateCountFile2(oneGenus,countFilePath)
    #cluster
    muthurParaStr="\"#cluster(column=%s, method=average, count=%s)\""%(distanceResult,countFilePath)
    commands.getoutput(muthurPath+" "+muthurParaStr)
    return True

def assignedSingleton2Annotated(seedDict,distanceFile,threshold,singletonfile):
    singletonList=[] #return 
    reDict={}
    usesingleton=[]
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
    usage="""usage: 
 --assigned/-a	(required) the file containing these taxonomically assigned sequences ("unique_sequences_multiple_assigned.fa")

 --unassigned/-u	(required) the file containing the sequences failing to be taxonomically assigned ("unique_sequences_multiple_unassigned.fa") 

 --taxonomy/-t	(required) the file of taxonomical information for all multiple sequences ("unique_sequences_multiple.taxonomy")

 --singleton/-s	(required) the file containing singleton tags ("unique_sequences_singleton.fa")

 --threshold/-r	(optional) to specify the distance cutoff for OTUs clustering, default:0.03

 --kdistance/-k	(optional) to specify cutoff of kmer distance, default:0.5

 --processors/-p	(optional) processors, default:1

 --help/-h	
    """
    processors=1
    threshold=0.03
    KmerFilter=0.5
    oldWorkDir=os.getcwd()
    opts,arg=getopt.getopt(sys.argv[1:],"a:u:t:s:r:p:k:h",['assigned=','unassigned=','taxonomy=','singleton=','threshold=','processors=','kdistance=','help',],)

    parameters=[a[0] for a in opts]
    if '-h' in parameters or '--help' in parameters:
        print usage
        sys.exit(1)
    if len(parameters)==0:
        print usage
        sys.exit(1)
    if '-a' not in parameters and '--assigned' not in parameters:
        print "***Error, --assigned/-a is requred.***\n"
        print usage
        sys.exit(1)
    if '-u' not in parameters and '--unassigned' not in parameters:
        print "***Error, --unassigned/-u is requred.***\n"
        print usage
        sys.exit(1)
    if '-t' not in parameters and '--taxonomy' not in parameters:
        print "***Error, --taxonomy/-t is requred.***\n"
        print usage
        sys.exit(1)
    if '-s' not in parameters and '--singleton' not in parameters:
        print "***Error, --singleton/-s is requred.***\n"
        print usage
        sys.exit(1)

    for i,a in opts:
        if i in ("--assigned","-a"):
            if not os.path.isfile(a):
                print "%s is not found."%(a)
                sys.exit(1)
            seedFasta=os.path.abspath(a)
        if i in ("--unassigned","-u"):
            if not os.path.isfile(a):
                print "%s is not found."%(a)
                sys.exit(1)
            sourceFasta=os.path.abspath(a)
        if i in ("--taxonomy","-t"):
            if not os.path.isfile(a):
                print "%s is not found."%(a)
                sys.exit(1)
            taxonomy=os.path.abspath(a)
        if i in ("--singleton","-s"):
            if not os.path.isfile(a):
                print "%s is not found."%(a)
                sys.exit(1)
            singletonfile=os.path.abspath(a)
        if i in ("--threshold","-r"):
            try:
                threshold=float(a)
            except:
                print "***Error, the value of threshold (--threshold/-r) must be decimals.***\n"
                print usage
                sys.exit(1)
        if i in ("--processors","-p"):
            try:
                processors=int(a)
            except:
                print "***Error, the processors (--processors/-p) must be integer.***\n"
                print usage
                sys.exit(1)
        if i in ("--kdistance","-k"):
            try:
                KmerFilter=float(a)
            except:
                print "***Error, the kmer distance (--kdistance/-k) must be decimals.***\n"
                print usage
                sys.exit(1)

    script_loc=os.path.split(os.path.realpath(sys.argv[0]))[0]
    remove_intermediate_file=r"rm -rf "+script_loc+r"/temp* "
    commands.getoutput(remove_intermediate_file)

    tempDir=script_loc+"/temp"+str(random.randint(10000,99999))
    os.mkdir(tempDir)


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


    candidateSeq=oldWorkDir+"/pengding_sequences_multiple.fa"
    singletonFileOut=oldWorkDir+"/pending_sequence_singleton.fa"
    outputOTUinfo=oldWorkDir+"/taxonomy_guided_OTU.list"
    outputOTUseq=oldWorkDir+"/taxonomy_guided_OTU.fa"
    outputTAX=oldWorkDir+"/taxonomy_guided_OTU.genus"

    #calculate distance of seed and candidate
    annotated_unannotated_distances=distanceCalculate(seedFasta,sourceFasta,KmerFilter,processors,tempDir,script_loc)
    seedDict,candidateList=clustSeedCandidate(seedFasta,annotated_unannotated_distances,threshold,sourceFasta)

    #assign singleton sequences to annotated sequences
    singleton_annotated_distances=distanceCalculate(seedFasta,singletonfile,KmerFilter,processors,tempDir,script_loc)
    seedDict,singletonList=assignedSingleton2Annotated(seedDict,singleton_annotated_distances,threshold,singletonfile)


    #calculate distance between seed
    SameGenusDict=createSeedTaxonomyDict(taxonomy,seedDict)
    seedSeqDict=fasta2dict(seedFasta)
    distanceEXE=script_loc+"/lib/distanceCalculateWithKmerOneFile"
    os.chdir(tempDir)
    k=0
    pool=multiprocessing.Pool(processes=processors)
    calculateGenus=[]
    bigGunus=[]
    for genus in SameGenusDict:
        if len(SameGenusDict[genus])<100:
            k+=1
            fastaFile=tempDir+"/OneGenus_%d.fa"%(k)
            pool.apply_async(multiOTUcall,(fastaFile,SameGenusDict[genus],distanceEXE,muthurPath))
        else:
            bigGunus.append(genus)
    pool.close()
    pool.join()
    seedseqDict=fasta2dict(seedFasta)
    for genus in bigGunus:
        k+=1
        multprocessForOneGenus(seedseqDict,SameGenusDict[genus],distanceEXE,muthurPath,processors,k,tempDir,script_loc)
    #conbine otu
    taoxnomyDict={}
    for line in open(taxonomy,"r"):
        lineList=line.strip().split()
        taoxnomyDict[lineList[0]]=lineList[1]


    totalSeqDict={}
    totalSeqDict.update(seedseqDict)
    totalSeqDict.update(fasta2dict(sourceFasta))
    totalSeqDict.update(fasta2dict(singletonfile))
    foutOTU=open(outputOTUinfo,'w')
    foutTax=open(outputTAX,'w')
    foutOTUseq=open(outputOTUseq,'w')
    OTUnum=1
    while k>0:
        fastaFile=tempDir+"/OneGenus_%d.fa"%(k)
        OTUlistFile=tempDir+"/OneGenus_%d.fa.an.unique_list.list"%(k)
        isouttax=False
        if os.path.isfile(OTUlistFile):
            OTUlistDict=mothurOTUlist2dict(OTUlistFile)
            thre=threshold
            while thre>=0:
                OTU=OTUlistDict.get(str(thre),None)
                if OTU:
                    OTUlist=OTU
                    for OTUone in OTUlist:
                        reOTU=[]
                        for seedName in OTUone.split(','):
                            if not isouttax:
                                otutaxnomy=taoxnomyDict[seedName]
                                isouttax=True
                            reOTU.append(seedName+"*")
                            candiList=seedDict.get(seedName,None)
                            reOTU+=candiList
                        OTUstr=",".join(reOTU)
                        foutOTU.write("%s%d\t%s\n"%("OTU",OTUnum,OTUstr))
                        foutTax.write("%s%d\t%s\n"%("OTU",OTUnum,otutaxnomy))
                        for OTUseqname in reOTU:
                            seqname=OTUseqname.strip('*')
                            foutOTUseq.write(">%s\n%s\n"%(seqname,totalSeqDict[seqname]))
                        OTUnum+=1
                    break
                else:
                    thre-=0.01
        else:
            for line in open(fastaFile,'r'):
                if line.startswith(">"):
                    seedName=line.strip().strip(">")
                    otutaxnomy=taoxnomyDict[seedName]
                    OTU=[seedName+"*"]
                    candiList=seedDict.get(seedName,None)
                    OTU+=candiList
            OTUstr=",".join(OTU)
            foutOTU.write("%s%d\t%s\n"%("OTU",OTUnum,OTUstr))
            foutTax.write("%s%d\t%s\n"%("OTU",OTUnum,otutaxnomy))
            OTUnum+=1
            for OTUseqname in OTU:
                seqname=OTUseqname.strip('*')
                foutOTUseq.write(">%s\n%s\n"%(seqname,totalSeqDict[seqname]))
        k-=1
    outputRetainedCandidate(sourceFasta,candidateList,candidateSeq)
    outputRetainedCandidate(singletonfile,singletonList,singletonFileOut)
    print "\noutput file:\n"
    print os.path.basename(outputOTUinfo)
    print os.path.basename(singletonFileOut)
    print os.path.basename(candidateSeq)
    print os.path.basename(outputTAX)

    remove_intermediate_file=r"rm -rf "+script_loc+r"/temp* "
    commands.getoutput(remove_intermediate_file)
