'''this script is used to convert Orthogroup tsv file from OrthoFinder2 into gene ID file per line'''
import re
import sys


def ParseOrthogroupTsv(filename):
    '''this function is used to parse the Orthogroup file into gene ID file,
    Note: gene ID file contains gene ID per ID a line'''

    geneids = []
    with open(filename) as input:
        for line in input:
            if 'Orthogroup' in line:
                continue
            else:
                parsed_line = re.split(', |\t', line.strip('\n'))
                # print(parsed_line)
                for x in parsed_line:
                    if not ('OG' in x or x == ''):
                        geneids.append(x)
    return geneids


geneids = ParseOrthogroupTsv(sys.argv[1])   # specify Orthogroup.tsv file
# print(geneids)

def LineOutput(geneids, outputfile):
    '''this function is used to rewrite the geneids'''
    output = open(outputfile, 'w')
    for id in geneids:
        output.write(id.replace('t', 'g') + '\n')
    output.close()


LineOutput(geneids, sys.argv[2])  # specify output filename