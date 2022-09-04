'''
Description: this script is used to parse the SNP extracted in a specific format, based on corresponding section
1) SNP format: "1",637325,chr08:55372944,55372944
We only need "chr08:55372944"
2) GFF file, 
e.g
LG01	EVM	gene	56614	64648	.	+	.	ID=Chr01G000010
Note: chromosome is named after "<LG><number>"
'''
from collections import defaultdict
import os
import sys


def ReadGFF(filename):
    '''this function is used to parse the GFF file'''
    info_list = defaultdict(list)
    with open(filename, 'rt') as input:
        for line in input:
            parsed_line = line.strip('\n').split('\t')

            name = parsed_line[0].replace('LG', 'chr')
            if 'Contig' in name: # only save the scaffold level chromosome
                continue
            # print(parsed_line[3])
            # print(parsed_line[4])
            info_list[name].append([int(parsed_line[3]), int(parsed_line[4]), parsed_line[len(parsed_line)-1].strip('ID=')])
    return info_list


def Readcsv(filename, info_list, number):
    '''this function is used to parse the specific SNP text file'''
    output = open('bio' + number + '.txt', 'w')  # specified output file name
    with open(filename, 'rt') as input:
        for line in input:
            name = line.split(',')[2].split(':')[0]
            pos = line.split(',')[2].split(':')[1]
            # print(type(pos))
            keys = list(info_list.keys())
            # print(keys)
            # print(name)
            for x in keys:
                if name == x:
                    for y in info_list[x]:
                        if int(pos) >= y[0] and int(pos) <= y[1]:
                            b_name = name.replace('chr', 'LG')
                            output.write(f'{b_name}\t{y[2]}\n')
                
    output.close()
                    

def main(biopath):
    '''
    main function to filter SNP in a run
    Note: should copy all input SNP file in the same dir
    '''
    info_list = ReadGFF('gene.txt')  # specify input GFF file
    for x in sorted(os.listdir(biopath), key = lambda i:i.split('bio')[1].split('.')[0]):
        # print(x.split('bio')[1].split('.')[0])
        Readcsv(os.path.join(biopath, x), info_list, x.split('bio')[1].split('.')[0])
        

main(sys.argv[1])