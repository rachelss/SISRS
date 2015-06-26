#!/usr/bin/env python2
from math import *
from collections import Counter
from Bio import SeqIO
import sys

#######################################
def print_matrix(matrix,seq1,seq2,weight=1):
    print '\t'+('\t'.join(map(str,list(seq1))))
    i = 0
    for line in matrix:
        print seq2[i]+"\t"+('\t'.join(map(str,line)))
        i +=1

def getmatrix(seq1,seq2,weight=1):
    indelValue = -1000
    matchValue = 500
    mismatchValue = indelValue
    
    matrix = []
    path = []
    rows = len(seq2) 
    cols = len(seq1)
    seq1 = "^"+seq1
    seq2 = "^"+seq2
    
    #make matrix and path map
    for i in range(rows+1):
        matrix.append([0]*(cols+1))
        path.append(["N"]*(cols+1))
        
    for i in range(1,rows+1):
        for j in range(1,cols+1):
            # penalty map
            if i==1 and j==1:
                from_left = matrix[i][j-1]
                from_top = matrix[i-1][j]
            elif i==rows and j==cols:
                from_left = matrix[i][j-1]
                from_top = matrix[i-1][j]
            elif i==1:
                from_left = matrix[i][j-1] + indelValue
                from_top = matrix[i-1][j]
            elif i==rows:
                from_left = matrix[i][j-1]
                from_top = matrix[i-1][j] + indelValue*weight   
            elif j==1:
                from_top= matrix[i-1][j] + indelValue*weight
                from_left = matrix[i][j-1]
            elif j==cols:
                from_top= matrix[i-1][j]
                from_left = matrix[i][j-1] + indelValue
             
            else:
                from_left = matrix[i][j-1] + indelValue
                from_top = matrix[i-1][j] + indelValue*weight
            if seq2[i]==seq1[j]:
                from_diag = matrix[i-1][j-1] + matchValue
            else:
                from_diag = matrix[i-1][j-1] + mismatchValue

            matrix[i][j]= max(from_left,from_top,from_diag)
            # path map
            if matrix[i][j]==from_left:
                path[i][j]="-"
            elif matrix[i][j]==from_top:
                path[i][j] = "|"
            elif matrix[i][j] == from_diag:
#                if 
                path[i][j] = "M"
            else:
                pass

#            if matrix[i][j]<0:
#                matrix[i][j]=0
    
    #print_matrix(matrix,seq1,seq2)
    #print 
    #print_matrix(path,seq1,seq2)
    
    return matrix,path

def getAlignment(seq1,seq2,matrix,path):
    addafter=[]
    seq1 = "^"+seq1
    seq2 = "^"+seq2
    iRow = len(matrix)-1
    jCol = len(matrix[0])-1
    #print iRow,jCol
    
    #trace correct path
    while iRow>=0:
        maxPnltyVaue = max(matrix[iRow])
            
        while jCol>=0:
            if matrix[iRow][jCol]==maxPnltyVaue:
                if iRow==len(matrix)-1 and jCol < len(matrix[0])-1:
                    addafter = seq1[jCol+1:]
                ipair = iRow
                jpair = jCol
                report2=[]
                report1=[]
                while (1):
                    #print "".join(report1)
                    #print "".join(report2)
                    if ipair ==0 and jpair==0:
                        break
                    if path[ipair][jpair]=="M":
                        report2.append(seq2[ipair])
                        report1.append(seq1[jpair])
                        ipair -= 1
                        jpair -= 1
                    elif path[ipair][jpair]=="-":
                        report2.append("-")
                        report1.append(seq1[jpair])
                        jpair -= 1
                    elif path[ipair][jpair]=="|":
                        report2.append(seq2[ipair])
                        report1.append("-")
                        ipair -= 1
                    elif path[ipair][jpair]=="N":
                        #print "N"
                        if ipair > 0:
                            report2.append(seq2[ipair])
                            report1.append("-")
                            ipair -=1
                        if jpair > 0:
                            report1.append(seq1[jpair])
                            report2.append("-")
                            jpair -= 1
                    #print 'row'+str(ipair)+' col'+str(jpair)+"\n"

            jCol -=1
        iRow -=1
        report1=report1[::-1]
        report2=report2[::-1]
        if addafter is not None:
            if len(addafter)==len(seq1)-1:
                report1=list(addafter)+report1
                report2=(['-']*len(addafter))+report2
            else:
                report1.extend(addafter)
                report2.extend(['-']*len(addafter))
        return report1,report2

def getMSA(s):
    if len(s)>1:
        #align first two reads
        matrix,path=getmatrix(s[0],s[1])
        a,b=getAlignment(s[0],s[1],matrix,path)
        s[0]="".join(a)
        s[1]="".join(b)
        pileup=[[]]*len(s[0])
       
        for i in range(len(s[0])):
            pileup[i]=[s[0][i]]     #position i sequence j
            pileup[i].append(s[1][i])
        i=2
        
        #align remaining reads
        while i<len(s):
            print i,
            cons = []
            for j in pileup:
                bases=[b for b in j if b in ['A','C','G','T','a','c','g','t']]
                if not bases:
                    print sys.argv[1], j
                    sys.exit(1)
                consbase=Counter(bases).most_common(1)[0][0]
                cons.append(consbase)
            matrix,path=getmatrix("".join(cons),s[i],len(pileup[0]))
            adjcons,b=getAlignment("".join(cons),s[i],matrix,path)
            s[i]="".join(b)
            #compare adjcons to cons and adjust pileup
            for j,b in enumerate(adjcons):
                if b=='-':
                    pileup.insert(j,['-']*i)
            for j in range(len(s[i])):    #j is position in sequence
                pileup[j].append(s[i][j])
            
            j+=1
            while j<len(pileup):
                pileup[j].append('-')
                j+=1
            i+=1
        
        allseqs=[]
        for i in range(len(pileup[0])):
            allseqs.append("".join([k[i] for j,k in enumerate(pileup) ]))
    else:
        return s[0]
    return allseqs

########################################

#get data
handle = open(sys.argv[1], "rU")
allseqs = list(SeqIO.parse(handle, "fasta"))
handle.close()

allseq2 = [str(s.seq) for s in allseqs ]    #just the sequence (as a string)
#lens = [len(s) for s in allseq2]            
#for i in range(max(lens)):
#    checkforN=[seq[i] for seq in allseq2 if i<len(seq)]
#    if len(set(checkforN))==1 and checkforN[0]=='N':
#        continue
#    else:
#        break
#for j,seq in enumerate(allseq2):
#    allseq2[j]=str(seq[i:])
#
#lens = [len(s) for s in allseq2]
#i=max(lens)-1
#while i>0:
#    checkforN=[seq[i] for seq in allseq2 if i<len(seq)]
#    if len(set(checkforN))==1 and checkforN[0]=='N':
#        i=i-1
#    else:
#        break
#        
#for j,seq in enumerate(allseq2):
#    if len(seq)>=i:
#        allseq2[j]=str(seq[:i])
#    elif len(allseq2[j])<i:
#        addto=['N']*(i-len(allseq2[j]))
#        allseq2[j]=list(allseq2[j])
#        allseq2[j].extend(addto)
#        allseq2[j]="".join(allseq2[j])

#lens = [len(s) for s in allseq2]
#for i in range(max(lens)):
#    checkforN=[seq[i] for seq in allseq2 if i<len(seq)]
#    if len(set(checkforN))==1 and checkforN[0]=='N':
#        print i

seqs=getMSA(allseq2)

folder_name='/'.join(sys.argv[1].split('/')[:-1])
filename=sys.argv[1].split('/')[-1]
outfile = open(folder_name+'/'+filename.replace('.fa','_align.fa'),'w')
for i,a in enumerate(allseqs):
    outfile.write('>'+a.name+"\n")
    outfile.write(seqs[i]+"\n")
outfile.close()
