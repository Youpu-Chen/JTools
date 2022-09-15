'''this module is used to extract gene id and calculate gene length'''
from Bio import SeqIO
import sys

def ParseFasta(filename):
    '''
    Note: the input file version should be like,
    A.gene2species.txt
    '''
    # Retrieve species name
    name = filename.split('.')[0]

    input = open(filename, 'r')
    
    info = []
    for record in SeqIO.parse(input, 'fasta'):
        id = record.id
        g_len = len(str(record.seq))
        info.append([id, str(g_len)])
        # 
        # print([id, g_len])
    
    for i in range(len(info)):
        if name in ['A', 'B', 'C']:
            print('\t'.join(info[i]) + '\t' + name + '.subgenome.w60')
        else:    
            print('\t'.join(info[i]) + '\t' + name)  # '\n' comes with print
    input.close()

ParseFasta(sys.argv[1])
