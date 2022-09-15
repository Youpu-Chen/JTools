'''this script is used to filter the ortholog based on the homologous chromosome info'''
from Bio import SeqIO
from collections import defaultdict
import re
import os
import sys


# class Ortholog():
#     def __init__(self, id):
#         self.id = id


def ParseInfo(filename):
    '''read the homologous info, and parse as dict'''
    pair_list = []
    with open(filename, 'rt') as input:
        for line in input:
            # print(line.strip('\n'))
            parsed_line = line.rstrip('\n').split('\t')
            
            tmp_list = []
            for x in parsed_line:
                if x == 'NA':
                    continue 
                else:
                    tmp_list.append(x)
            
            pair_list.append(tmp_list)
    pair_list = sorted(pair_list, key = lambda x:int(x[0][1:]))  # using chromosome number of first subgenome to sort the list
    # for x in pair_list:
    #     print(x)

    return pair_list
            

def ParseOrthologFasta(filename):
    '''this function is used to parse the Ortholog fasta from OrthoFinder2
    Note: OrthoFinder2 will read the faa in the working dir, and sort the filename automatically
    My Design (2022-09-06 version): 
    '''
    input = open(filename,'rt')
    pattern = re.compile(r'[A-C][\d]+.')
    # ortho_dict = defaultdict(list)
    ortho_dict = {}
    id_list = []
    seq_list = []
    for record in SeqIO.parse(input, 'fasta'):

        if re.search(pattern, record.id):  # only save sequence from the A, B, C subgenome
            id_list.append(record.id)
            seq_list.append(str(record.seq))
    
    # sort
    ortho_dict[','.join(id_list)] = seq_list

    return ortho_dict


def FilterOrtholog(pair_list, ortho_dict, outputfile):
    '''this function is used to filter the Ortholog fatsa'''
    # 
    output = open(outputfile, 'w')
    # 
    pattern = re.compile(r'[A-C][\d+]')
    
    # Reform the keys
    # print(type(ortho_dict.keys()))
    matches = re.findall(pattern, list(ortho_dict.keys())[0])
    print(f'the re matches is {matches}')

    # Filter out all the uncomplete or uncorrect-paired orthologs
    split_ortho_keys = list(ortho_dict.keys())[0].split(',')
    print(split_ortho_keys)
    ortho_values = list(ortho_dict.values())
    # print(f'the ortho_values is {ortho_values}')
    # print(f'the type of ortho_values is {type(ortho_values)}')

    # Rebuild the pair list
    tmp_list = []
    print(f'the pair_list is {pair_list}')
    for x in pair_list:
        tmp_list.append(','.join(x).rstrip(','))
    print(f'The tmp list is {tmp_list}')

    if ','.join(matches) in tmp_list:
        print(','.join(matches))
        for i in range(len(ortho_values[0])):
            print(i)
            output.write(split_ortho_keys[i] + '\n')
            output.write(ortho_values[0][i] + '\n')
    output.close()
            
        
def main(dir):
    input_list = os.listdir(dir)
    pair_list = ParseInfo(sys.argv[1])   # specify the chr_pair.txt

    for fa in input_list:
        ortho_dict = ParseOrthologFasta(os.path.join(dir, fa))
        FilterOrtholog(pair_list, ortho_dict, os.path.join(dir, fa.split('.')[0] + '.filtered.fa'))

main(sys.argv[2])  # specify the Ortholog fasta directory


if __name__ == "__main__":
    # Test
    # pair_list = ParseInfo('../info/chr_pair.txt')
    # # print(pair_list)
    # ortho_dict = ParseOrthologFasta('../test_data/OG0005241.fa')
    # # print(ortho_dict)

    # FilterOrtholog(pair_list, ortho_dict, '../test_data/filter.fa')


    # 
    pass
