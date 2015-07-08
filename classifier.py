from subprocess import call
import os,sys,getopt
import commands
import linecache
import stat
from sys import path
path.append(sys.path[0]+"/lib/")
import FuncbioOTU
from FuncbioOTU import fasta2dict

def getSeqName(taxonomyFile):
    levelValue=-(0+2)
    seedSeqNameList=[]
    candidateSeqNameList=[]
    for line in open(taxonomyFile,'r'):
        SequenceName=line.split()[0]
        line_list=line.strip(";").split(";")
        targetLevel=line_list[levelValue]
        if targetLevel=="unclassified":
            candidateSeqNameList.append(SequenceName)
        else:
            seedSeqNameList.append(SequenceName)
    return seedSeqNameList,candidateSeqNameList

def getSeedSequence(seedSeqNameList,sourceFasta,seedFasta):
    fseed=open(seedFasta,"w")
    sourceDict=fasta2dict(sourceFasta)
    for sequenceName in seedSeqNameList:
        seedSequence=sourceDict.get(sequenceName,False)
        if seedSequence:
            fseed.write(">%s\n"%(sequenceName))
            fseed.write("%s\n"%(seedSequence))
        else:
            print "Error, % not found."
            sys.exit(1)
    fseed.close()
    return True

def getCandidateSequence(candidateSeqNameList,sourceFasta,candidateFasta):
    fcand=open(candidateFasta,"w")
    sourceDict=fasta2dict(sourceFasta)
    for sequenceName in candidateSeqNameList:
        candSequence=sourceDict.get(sequenceName,False)
        if candSequence:
            fcand.write(">%s\n"%(sequenceName))
            fcand.write("%s\n"%(candSequence))
        else:
            print "Error, % not found."
            sys.exit(1)
    fcand.close()
    return True

def getSeedMain(taxonomyFile,sourceFasta,seedFasta,candidateFasta):
    seedSeqNameList,candidateSeqNameList=getSeqName(taxonomyFile)
    getSeed=getSeedSequence(seedSeqNameList,sourceFasta,seedFasta)
    getCandidate=getCandidateSequence(candidateSeqNameList,sourceFasta,candidateFasta)
    return True
if __name__=="__main__":

    usage="""usage:
    --fasta/-f   (required) the input file (in fasta format) that contains all your sequences for taxonomical assignments ("unique_sequences_multiple.fa")

    --template/-t   (required) the reference database, such as that are retrieved from RDP or SILVA (in fasta format)

    --taxonomy/-a   (required) the taxonomical file corresponding to your reference sequences

    --ksize/-k   (optional) to specify the size of kmer, default:8

    --iters/-i   (optional) to specify the number of iteration, default:1000

    --cutoff/-c   (optional) to specify bootstrap cutoff for supporting the taxonomical assignment, default:95

    --processors/-p   (optional) processors, default:1
 """

    opts,arg=getopt.getopt(sys.argv[1:],"f:t:a:m:k:i:c:n:s:h:b:g:x:p:",['fasta=','template=','taxonomy=','method=','ksize=','iters=','cutoff=','numwanted=','search=','match=','mismath=','gapopen=','gapextend=','processors='],)
#Required parameter check
    parameters=[a[0] for a in opts]
    paraDict={}
    paraDict["method"]="wang"
    paraDict["iters"]=1000
    paraDict["cutoff"]=95
    paraDict["processors"]=1
    if "--fasta" not in parameters and "-f" not in parameters:
        print usage
        sys.exit(1)        
    if "--template" not in parameters and "-t" not in parameters:
        print usage
        sys.exit(1)
    if "--taxonomy" not in parameters and "-a" not in parameters:
        print usage
        sys.exit(1)
    #optional parameter check
    for i,a in opts:
        if i in ("-f","--fasta"):
            if not os.path.isfile(a):
                print "***Error, %s is not found.***"%(a)
                sys.exit(1)
            fasta=a
            paraDict["fasta"]=a
        if i in ("-t","--template"):
            template=a
            paraDict["template"]=a
            if not os.path.isfile(a):
                print "***Error, %s is not found.***"%(a)
                sys.exit(1)
        if i in ("-a","--taxonomy"):
            if not os.path.isfile(a):
                print "***Error, %s is not found.***"%(a)
                sys.exit(1)
            paraDict["taxonomy"]=a

        if i in ("-k","--ksize"):
            try:
                paraDict["ksize"]=int(a)
            except:
                print "***Error, the value of ksize (--ksize/-k) must be integer.***\n"
                print usage
                sys.exit(1)
        if i in ("-i","--iters"):
            try:
                paraDict["iters"]=int(a)
            except:
                print "***Error, the value of iters (--iters/-i) must be integer.***\n"
                print usage
                sys.exit(1)
        if i in ("-c","--cutoff"):
            try:
                paraDict["cutoff"]=int(a)
            except:
                print "***Error, the value of cutoff (--cutoff/-c) must be integer.***\n"
                print usage
                sys.exit(1)
        if i in ("-p","--processors"):
            try:
                paraDict["processors"]=int(a)
            except:
                print "***Error, the processors (--processors/-p) must be integer.***\n"
                print usage
                sys.exit(1)
    #set probs to F
    paraDict["probs"]="F"
    useParaList=[]
    #optional
    for key in paraDict:
        value=paraDict[key]
        para="%s=%s"%(key,value)
        useParaList.append(para)
    #
    muthurParaStr="\"#classify.seqs("+",".join(useParaList)+")\""
    script_loc=os.path.split(os.path.realpath(sys.argv[0]))[0]
    muthurPath=script_loc+"/lib/Mothur.cen_64/mothur/mothur"
    uchime=script_loc+"/lib/Mothur.cen_64/mothur/uchime"
    try:
        os.chmod(muthurPath,stat.S_IRWXU)
    except:
        commands="chmod a+x %s"%(muthurPath)
        print "\n***Please give executable permission to %s by this terminal commands:\n%s  ***"%(muthurPath,commands)
        sys.exit(1)
    try:
        os.chmod(uchime,stat.S_IRWXU)
    except:
        commands="chmod a+x %s"%(uchime)
        print "\n***Please give executable permission to %s by this terminal commands:\n%s  ***"%(uchime,commands)
        sys.exit(1)
    returnMsg=commands.getoutput(muthurPath+" "+muthurParaStr)
    msyList=returnMsg.split("\n")
    msyList.reverse()
    for lastLine in msyList:
        if lastLine.endswith(".taxonomy") and os.path.isfile(lastLine):
            taxonomyPath=lastLine
            break
    for lastLine in msyList:
        if lastLine.endswith(".tax.summary") and os.path.isfile(lastLine):
            summaryPath=lastLine
            break

    workDir=os.path.dirname(os.path.abspath(fasta))
    basename=os.path.basename(fasta).split('.')[0]
    seedFastaPath=workDir+"/unique_sequences_multiple_assigned.fa"
    candidateFastaPath=workDir+"/unique_sequences_multiple_unassigned.fa"
    getSeedMain(taxonomyPath,fasta,seedFastaPath,candidateFastaPath)
    print "\noutput file:\n"
    print os.path.basename(seedFastaPath)
    print os.path.basename(candidateFastaPath)
    rmsummary="rm %s"%(summaryPath)
    commands.getoutput(rmsummary)
    taxonomyNewPath=workDir+"/unique_sequences_multiple.taxonomy"
    mvcommands="mv %s %s"%(taxonomyPath,taxonomyNewPath)
    commands.getoutput(mvcommands)
    print os.path.basename(taxonomyNewPath)
