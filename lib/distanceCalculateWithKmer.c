#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define kmerSize 7
#define hashSize 5*5*5*5*5*5*5

#define MISMATCH -1
#define MATCH 1
#define GAPOPEN -2
#define GAPEXTEND -1
struct kmerpoint{
    short isuse;
    short kmerA;
    short kmerB;
};

struct point{
    int score;  /*record the score of the point.*/
    int direction; /*record the direction where is the score from.0-diagonal, 1-top, 2-left, */
};

struct fasta{
    char name[200];
    char sequence[1000];
};

struct distance{
    char seed[200];
    char candidate[200];
    float dis;
};

int Hash(register char dna[]){
    register int sum;
    sum=0;
    register int digital;
    register int i;
    register int dnaLength=strlen(dna);
    register int n=1;
    for (i=0;i<dnaLength;i++){
       digital=((int)((int)dna[i]%10*0.5+0.5))*n;
       n*=5;
       sum+=digital;
    }
    return sum;
}

float kmerDistance(register char seq1[], register char seq2[]){
    struct kmerpoint hashList[hashSize];
    register int seq1Length=strlen(seq1);
    register int seq2Length=strlen(seq2);
    int hashValueL[2000];
    register int i,j,k=0;
    char kmer[kmerSize];
    register int a=0,b=0;
    register int key;
    register float distance=0;

    for(i=0;i<seq1Length-kmerSize+1;i++){
        for (j=0;j<kmerSize;j++){
            kmer[j]=seq1[i+j];
        }
        key=Hash(kmer);
        if (hashList[key].isuse==1){
            hashList[key].kmerA+=1;
            hashValueL[k]=key;
            k++;
        }
        else{
            hashList[key].isuse=1;
            hashList[key].kmerA=1;
            hashList[key].kmerB=0;
            hashValueL[k]=key;
            k++;
        }
    }

    for(i=0;i<seq2Length-kmerSize+1;i++){
        for (j=0;j<kmerSize;j++){
            kmer[j]=seq2[i+j];
        }
        key=Hash(kmer);
        if (hashList[key].isuse==1){
            hashList[key].kmerB+=1;
        }
        else{
            hashList[key].kmerB=1;
            hashList[key].kmerA=0;
            hashList[key].isuse=1;
            hashValueL[k]=key;
            k++;
        }
    }
   int minLen;
   int minkmer=0,sumkmer=0;
   minLen=seq1Length<seq2Length?seq1Length:seq2Length;
   for (i=0;i<k;i++){
        key=hashValueL[i];
        b=hashList[key].kmerB;
        a=hashList[key].kmerA;
        minkmer=a<b?a:b;
        sumkmer+=minkmer;
        hashList[key].isuse=0;
        hashList[key].kmerA=0;
        hashList[key].kmerB=0;
    }
    distance=1-((float)sumkmer/(minLen-kmerSize+1));
    return distance;
}


float Init(char seq1[], char seq2[]){  /*input a pointer for Init function*/
    struct point table[1000][1000]; 
    register int i,j,m,n;
    register int a_score,b_score,c_score;
    register int mydirect;
    register int myscore;
    register int penalty;
    register int lineLen=strlen(seq1)+2;
    register int rowLen=strlen(seq2)+2;
    char myseq1[1000];
    char myseq2[1000];
    register int NID=0;
    register float distance;
    register int mismatch=0,match=0,gap=0;
    /*initializing the first and second line*/
    for (i=0;i<lineLen;i++){
        if (i<2){
            table[0][i].score=0;
            table[0][i].direction=0;
            table[1][i].score=0;
            table[1][i].direction=0;
        }
        else{
            table[0][i].score=(int)seq1[i-2];
            table[0][i].direction=0;
            penalty=GAPOPEN+(i-2)*GAPEXTEND;
            table[1][i].score=penalty;
            table[1][i].direction=2;
        }
    }
    /*initializing the first and second row*/
    for (j=2;j<rowLen;j++){
        table[j][0].score=(int)seq2[j-2];
        table[j][0].direction=0;
        penalty=GAPOPEN+(j-2)*GAPEXTEND;
        table[j][1].score=penalty;
        table[j][1].direction=1;
    }

    /*fill table start*/
    for (j=2;j<rowLen;j++){
        for (i=2;i<lineLen;i++){
            /*check current point score and direction*/
            a_score=table[j-1][i-1].score;     /*a-dialogal, b-top, c--left*/
            if (table[j][0].score==table[0][i].score){

                myscore=a_score+MATCH;
                mydirect=0;
            }
            else{
                    myscore=a_score+MISMATCH;
                    mydirect=0;
            }

            /*check the top*/
            if (table[j-1][i].direction==1){
                b_score=table[j-1][i].score+GAPEXTEND;
                if (b_score>=myscore){
                    myscore=b_score;
                    mydirect=1;
                }
            }
            else{
                b_score=table[j-1][i].score+GAPOPEN;
                if (b_score>=myscore){
                    myscore=b_score;
                    mydirect=1;
                }
            }
            /*check the top*/
            if (table[j][i-1].direction==2){
                c_score=table[j][i-1].score+GAPEXTEND;
                if (c_score>=myscore){
                    myscore=c_score;
                    mydirect=2;
                }
            }
            else{
                c_score=table[j][i-1].score+GAPOPEN;
                if (c_score>=myscore){
                    myscore=c_score;
                    mydirect=2;
                }
            }
            /*fill the point information*/
            table[j][i].score=myscore;
            table[j][i].direction=mydirect;
        }
    }
 

    /*get alignment sequence start*/
    i=lineLen;
    j=rowLen;
    while (j>=1 && i>=1){
            mydirect=table[j][i].direction;
            if (mydirect==0){
                myseq1[NID]=table[0][i].score;
                myseq2[NID]=table[j][0].score;
                NID++;
                i--;
                j--;
            }
            else if (mydirect==1){
                myseq1[NID]='-';
                myseq2[NID]=(char)(table[j][0].score);
                j--;
                NID++;
            }
            else{
                myseq1[NID]=(char)(table[0][i].score);
                myseq2[NID]='-';
                i--;
                NID++;
            }

    }
    NID--;
    NID--;

    /*calculate distance start*/
    /*delete the start and end gap.*/
    j=NID;
    while (1){
        if (myseq1[j]=='-' || myseq2[j]=='-') j--;
        else break;
    }
    i=1;
    while (1){
        if (myseq1[i]=='-' || myseq2[i]=='-') i++;
        else break;  
    }
    /*calculate distance start*/
    while(j>=i){
        if (myseq1[j] == myseq2[j]){
            match++;
        }
        else if (myseq1[j] != '-' && myseq2[j] != '-'){
            mismatch++;
        }
        else if (myseq1[j-1] != '-' && myseq2[j-1] != '-'){
            gap++;
        }
        j--;
    }

    distance=(float)(gap+mismatch)/(float)(gap+mismatch+match);
    /*calculate distance end*/
    return distance;
}
 
void ReadFasta(char candidatefilepath[], char seedfilepath[],char output[],float kthreshould){
    int len;
    int cfaLen=0;
    int sfaLen=0;
    char *pc,*ps;
    int i,j;
    int seqID=-1;
    int a=0,b=0,c=0;
    int id=0;
    int csequence[1000];
    int cname[200];
    int ssequence[1000];
    int sname[200];
    struct fasta * candidate;
    struct fasta * seed;
    FILE *fa;
    fa=fopen(candidatefilepath,"rb");
    fseek(fa,0L, SEEK_END);
    len=ftell(fa);
    fseek(fa,0L,SEEK_SET);
    pc=(char*)malloc(len+1);
    memset(pc,0x00,len+1); 
    fread(pc,len,sizeof(char),fa);
    for (i=0;i<len;i++){
        if (pc[i]=='>'){
            cfaLen++;
        }
    }
    /*create a list for candidate fasta*/
    candidate = (struct fasta *)malloc(cfaLen*sizeof(struct fasta));
    for (i=0;i<len;i++){
        if (pc[i]=='>'){
            seqID++;
            a=1;
        }
        else if (pc[i] != '\n' && a==1){
            candidate[seqID].name[b]=pc[i];
            b++;
        }
        else if (pc[i] == '\n' && a==1){
            a=0;
            b=0;
        }

        else if(a==0 && pc[i] != '\n'){
            candidate[seqID].sequence[c]=pc[i];
            c++;
        }
        else if (a==0 && pc[i] == '\n'){
            c=0;            
        }
    }

    /*-----------for seed---------------*/
    fa=fopen(seedfilepath,"rb");
    fseek(fa,0L, SEEK_END);
    len=ftell(fa);
    fseek(fa,0L,SEEK_SET);
    ps=(char*)malloc(len+1);
    memset(ps,0x00,len+1); 
    fread(ps,len,sizeof(char),fa);
    for (i=0;i<len;i++){
        if (ps[i]=='>'){
            sfaLen++;
        }
    }
    /*create a list for seed fasta*/
    seqID=-1;
    a=b=c=0;
    seed = (struct fasta *)malloc(sfaLen*sizeof(struct fasta));
    for (i=0;i<len;i++){
        if (ps[i]=='>'){
            seqID++;
            a=1;
        }
        else if (ps[i] != '\n' && a==1){
            seed[seqID].name[b]=ps[i];
            b++;
        }
        else if (ps[i] == '\n' && a==1){
            a=0;
            b=0;
        }

        else if(a==0 && ps[i] != '\n'){
            seed[seqID].sequence[c]=ps[i];
            c++;
        }
        else if (a==0 && ps[i] == '\n'){
            c=0;            
        }
    }
    //printf("name:%s,sequence:%s",seed[sfaLen-1].name,seed[sfaLen-1].sequence);

   //calculate the distance of seed sequence and candidate sequence.
    register float rdistance=0;
    float kmerDis=0;
    FILE *fout;
    fout=fopen(output,"wb");
    //struct point table[1000][1000];
    for (i=0;i<sfaLen;i++){
        if (strlen(seed[i].sequence)<=1) break;
        for (j=0;j<cfaLen;j++){
            if (strlen(candidate[j].sequence)<=1) break;

            //calculate kmer distance at first
            kmerDis=kmerDistance(seed[i].sequence,candidate[j].sequence);
            if (kmerDis<=kthreshould){
                rdistance=Init(seed[i].sequence,candidate[j].sequence);
                fprintf(fout,"%s\t%s\t%f\t%f\n",seed[i].name,candidate[j].name,rdistance,kmerDis);
            }
        }

    }
    free(ps);
    free(pc);
    fclose(fout);
}

void main(int argc, char *argv[]){
    "usage: xxx seed/FILE/PAATH candidate/file/path output/file/path kthreshould";
    ReadFasta(argv[2],argv[1],argv[3],atof(argv[4]));
}
