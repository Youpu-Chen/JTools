'''
Description: 
- this file is used to parse the GFF file

Usage:
FastafromGFF.py input.gff input.fa outputfilename
Note: input.fa should only contain ID, not full header
e.g. 
>ID
not like, >ID xxxx xxxx
'''
from collections import defaultdict
from Bio import SeqIO
import sys


class GFFline:
    def __int__(self, chr, start, end, id):
        self.chr = chr
        self.start = start
        self.end = end
        self.id = id


class Fasta:
    def __init__(self, id, seq):
        self.id = id
        self.seq = seq


# def ParseAttributes(attributes):
#     '''this function is used to parse the GFF attributes'''
#     # Preset variables
#     sep = ';'
#     attributes_dict = {}
    
#     # 
#     parsed_attr = attributes.split(sep)
#     for x in parsed_attr:
#         tmp = x.split('=')
#         key = tmp[0]
#         value = tmp[1]
#         attributes_dict[key] = value
    
#     return attributes_dict


def ParseGFF(file):
    '''this function is used to read gff file and extract the start and end information'''
    # Preset
    # chr_attrs = defaultdict(list)
    line_list = defaultdict(list)
    
    # Parse the GFF
    input = open(file, 'rt')

    for line in input:
        if '#' in line:  # skip the comment infor
            continue
        else:
            parsed_line = line.strip('\n').split('\t')
            gffline = GFFline()
            gffline.chr = parsed_line[0]
            gffline.start = parsed_line[3]
            gffline.end = parsed_line[4]
            # chr_attrs[parsed_line[0]].append(ParseAttributes(parsed_line(len(parsed_line - 1))))
            print(gffline.chr)
            line_list[gffline.chr].append( [int(gffline.start), int(gffline.end)] )
    
    input.close()
    return line_list
            

def ParseFasta(file):
    '''this function is used to parse the genome fasta'''
    # Preset
    fasta_list = {}
    
    # Create instance
    with open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            fasta = Fasta(record.id, record.seq)
            # fasta.id = record.id
            # fasta.seq = record.seq
            fasta_list[fasta.id] = fasta.seq
    return fasta_list


def Retrieveseq(info, fasta, filename):
    '''this function is used to retrieve sequences from deposited GFF list
    Note:
    1) info -> line_list
    2) fasta -> fasta_list
    '''
    # 
    output = open(filename, 'w')

    # Build items and index
    c = sorted(fasta.items(), key=lambda x:x[0])
    index = [x[0] for x in c]

    # Loop through gff list to retrieve sequences
    for line in info.items():
        chr = line[0]
        i = index.index(chr)

        for x in line[1]:
            seq_name = f'{chr}_{str(x[0])}-{str(x[1])}'
            output.write('>' + seq_name + '\n')
            output.write(str( c[i][1][(int(x[0])-1):(int(x[1])-1)] ) + '\n')
            print('Seq length: %s' % ( (int(x[1])-1) - (int(x[0])-1) ))
    output.close()



def main():
    line_list = ParseGFF(sys.argv[1])
    fasta_list = ParseFasta(sys.argv[2])

    Retrieveseq(line_list, fasta_list, sys.argv[3])
    print('Seq retrieve completed!')


# def main():
#     line_list = ParseGFF('Arabidopsis_suecica_CDS.gff')
#     fasta_list = ParseFasta('test.fa')

#     Retrieveseq(line_list, fasta_list, 'CDS.fa')
#     print('Seq retrieve completed!')

if __name__ == "__main__":
    main()