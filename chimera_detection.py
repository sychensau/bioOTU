import sys
import os
import stat
import getopt
import random
import commands
import shutil
from sys import path
path.append(sys.path[0]+"/lib/")
from FuncbioOTU import getAttriValueFromSeqName,fasta2dict
def selectRefBySampeSize(minsampleSize,tempDir,copysequence,reference):
    totalRef=os.path.splitext(copysequence)[0]+".ref"
    print totalRef
    fref=open(totalRef,'w')  #used as reference for chimera detection
    requry=os.path.splitext(copysequence)[0]+".requry"
    fqury=open(requry,'w') #used as requery for chimera detecthion
    nonchimera=os.path.splitext(copysequence)[0]+".nonchimera"
    fnonch=open(nonchimera,'w') #a part of sequence for output that nonchimera.
    requrySeqDict=fasta2dict(copysequence)
    for name in requrySeqDict:
        if int(getAttriValueFromSeqName(name,"sampleSize"))>=minsampleSize:
            fref.write(">%s\n%s\n"%(name,requrySeqDict[name]))
            fnonch.write(">%s\n%s\n"%(name,requrySeqDict[name]))
        else:
            fqury.write(">%s\n%s\n"%(name,requrySeqDict[name]))
    refdict=fasta2dict(reference)
    for name in refdict:
        fref.write(">%s\n%s\n"%(name,refdict[name]))
    return totalRef,requry,nonchimera
if __name__=="__main__":
    usage="""usage: 

 --best_reference/-b	(required) input file containing all reference sequences for chimera detection ("taxonomy_guided_OTU.fa")

 --pengding_sequence/-s	(required) input file containing all candidate sequences subject to chimera detection ("pengding_sequences_multiple.fa")

 --sample_size/-r	(optional) to specify minimum value of sample size for adding these pengding sequences into reference database, default:1

 --processors/-p	(optional) processors, default:4
    """
    oldWorkDir=os.getcwd()
    useParaList=[]
    size='a'
    minsamplesize=False
    processors=4
    opts,arg=getopt.getopt(sys.argv[1:],"b:s:r:p:h",['best_reference=','pengding_sequence=','sample_size=','processors=','help'],)

    parameters=[a[0] for a in opts]
    if '-h' in parameters or '--help' in parameters:
        print usage
        sys.exit(1)
    if len(parameters)==0:
        print usage
        sys.exit(1)
    if '-b' not in parameters and '--best_reference' not in parameters:
        print "***Error, --best_reference/-b is requred.***\n"
        print usage
        sys.exit(1)
    if '-s' not in parameters and '--pengding_sequence' not in parameters:
        print "***Error, --pengding_sequence/-s is requred.***\n"
        print usage
        sys.exit(1)

    for i,a in opts:
        if i in ("--best_reference","-b"):
            if not os.path.isfile(a):
                print "%s is not found."%(a)
                sys.exit(1)
            reference=os.path.abspath(a)
        if i in ("--pengding_sequence","-s"):
            if not os.path.isfile(a):
                print "%s is not found."%(a)
                sys.exit(1)
            sequence=os.path.abspath(a)
        if i in ("--sample_size","-r"):
            try:
                minsamplesize=int(a)
            except:
                print "***Error, sample size (--sample_size/-r) must be integer.***\n"
                print usage
                sys.exit(1)
        if i in ("--processors","-p"):
            try:
                processors=int(a)
            except:
                print "***Error, the processors (--processors/-p) must be integer.***\n"
                print usage
                sys.exit(1)

    nonchimera=oldWorkDir+"/pengding_sequences_multiple.nonchimera"
    chimera=oldWorkDir+"/pengding_sequences_multiple.chimera"
    script_loc=os.path.split(os.path.realpath(sys.argv[0]))[0]
    tempDir=script_loc+"/temp"+str(random.randint(10000,99999))
    os.mkdir(tempDir)
    os.chdir(tempDir) #change work directory

    copysequence=tempDir+"/"+os.path.basename(sequence)
    shutil.copyfile(sequence,copysequence)

    uchimeraPath=script_loc+"/lib/Mothur.cen_64/mothur/uchime"
    mothurPath=script_loc+"/lib/Mothur.cen_64/mothur/mothur"
    try:
        os.chmod(uchimeraPath,stat.S_IRWXU)
    except:
        commands="chmod a+x %s"%(uchimeraPath)
        print "Please give executable permission to %s by this terminal commands:\n%s"%(uchimeraPath,commands)
    try:
        os.chmod(mothurPath,stat.S_IRWXU)
    except:
        commands="chmod a+x %s"%(mothurPath)
        print "Please give executable permission to %s by this terminal commands:\n%s"%(mothurPath,commands)

        sys.exit(1)
    if minsamplesize:
        totalRef,requry,subnonchimera=selectRefBySampeSize(minsamplesize,tempDir,copysequence,reference)
        useParaList.append("fasta=%s"%(requry))
        useParaList.append("reference=%s"%(totalRef))
    else:
        useParaList.append("fasta=%s"%(copysequence))
        useParaList.append("reference=%s"%(reference))

    useParaList.append("processors=%s"%(processors))

    muthurParaStr="\"#chimera.uchime("+",".join(useParaList)+")\""
    muthurPath=script_loc+"/lib/Mothur.cen_64/mothur/mothur"
    commands.getoutput(muthurPath+" "+muthurParaStr)

    #put chimera sequence names into list. 
    accons=os.path.splitext(copysequence)[0]+".uchime.accnos"
    acconsList=[]
    for line in open(accons,"r"):
        acconsList.append(line.strip())

    faDict=fasta2dict(copysequence)

    fout1=open(nonchimera,"w")
    fout2=open(chimera,"w")

    for line in faDict:
        if line in acconsList:
            fout2.write(">%s\n%s\n"%(line,faDict[line]))
        else:
            fout1.write(">%s\n%s\n"%(line,faDict[line]))
    if minsamplesize:
        subnonDict=fasta2dict(subnonchimera)
        for name in subnonDict:
            fout1.write(">%s\n%s\n"%(name,subnonDict[name]))

    remove_intermediate_file=r"rm -rf "+script_loc+r"/temp* "
    commands.getoutput(remove_intermediate_file)
