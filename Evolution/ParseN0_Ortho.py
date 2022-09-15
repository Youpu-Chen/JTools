'''this module is used to parse the N0.tsv from OrthoFinder2 to corresponding list'''
import re
import sys


def ParseN0(filename, outputfile):
    '''this function is used to parse the N0.tsv'''
    output = open(outputfile, 'w')
    with open(filename, 'rt') as input:
        next(input)
        for line in input:
            parsed_line = re.split(', |\t', line.strip('\n'))
            # remove empty string 
            while "" in parsed_line:
                parsed_line.remove("")
            print(parsed_line)
            for x in range(3, len(parsed_line)):
                output.write(f"{parsed_line[0]}\t{parsed_line[x].replace('t', 'g')}\n")


ParseN0(sys.argv[1], sys.argv[2])   # specify N0.tsv and output filename