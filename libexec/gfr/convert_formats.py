#! /usr/local/bin/python
import sys
from Bio import AlignIO,SeqIO
from Bio.Alphabet import IUPAC

def wrap(text, width):
    lines=len(text)/width
    lines+=1
    newtext=[list() for i in range(lines)]
    z=0
    for i,j in enumerate(text):
        if len(newtext[z])<width:
            newtext[z]+=j
        else:
            z+=1
            newtext[z]+=j
    for i,j in enumerate(newtext):
        newtext[i]="".join(j)
    return newtext
    
    
####################################
def main():
    if len (sys.argv) != 4 :
        print "Please provide file, the file format, and the desired file format "
        sys.exit (1)
    else:
        f = sys.argv[1]
        fout = "".join(f.split('.')[:-1])
        formatin = sys.argv[2]
        formatout  = sys.argv[3]
        if formatout == 'nexus':
            AlignIO.convert(f,formatin,fout+'.'+formatout,formatout,alphabet= IUPAC.ambiguous_dna)
        if formatout == 'mega':
            handle = open(f, "rU")
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "phylip-relaxed"))
            handle.close()
            
            outfile = open(fout+'.'+formatout,'w')
            outfile.write('#mega'+"\n")
            outfile.write('!Title Mytitle;'+"\n")
            outfile.write('!Format DataType=DNA indel=-;'+"\n\n")
            
            for n in record_dict:
                outfile.write('#'+n+"\n")
                newseq=wrap(str(record_dict[n].seq),60)
                for s in newseq:
                    outfile.write(s+"\n")
            
            outfile.close()
        else:
            AlignIO.convert(f,formatin,fout+'.'+formatout,formatout)
if __name__ == "__main__":
    main() 