'''this script is used to seek the gene belonging relationships
Note:
1) Orthogroup.tsv
2) Blast input
'''
from collections import defaultdict
import re
import os
import sys

def ParseOrthogroup(filename):
    
    orthos = defaultdict(list)

    with open(filename, 'rt') as input:
        next(input)
        for line in input:
            parsed_line = re.split(', |\t', line.strip())
            while '' in parsed_line:
                parsed_line.remove('')
            
            orthos[parsed_line[0]] = parsed_line[1:]
    return orthos


def ParseBlast(filename):
    with open(filename, 'rt') as input:
        line = input.readline()
        query_id = line.split('\t')[0]
        ref_id = line.split('\t')[1]
    return query_id, ref_id


def main():
    dir = sys.argv[1]  # specify the blast input path
    ref_ids = []
    for i in os.listdir(dir):
        query_id, ref_id = ParseBlast(os.path.join(dir, i))
        ref_ids.append(ref_id)
    # print(ref_ids)
    # print(len(ref_ids))
    joined_ref_ids = '\t'.join(ref_ids)

    orthos = ParseOrthogroup(sys.argv[2]) # specify Orthologgroup.tsv
    ortho_number, orthos_context  = [x[0] for x in orthos.items()], [x[1] for x in orthos.items()]
    # Loop through the items to get the Orthogroup votes
    votes = {}
    for ref_id in ref_ids:
        for index, context in enumerate(orthos_context):
            for ortho in context:
                if ref_id == ortho:
                    votes[ortho_number[index]] = votes.get(ortho_number[index], 0) + 1
    # print(votes)
    print(f'{query_id}\t{max(votes, key=votes.get)}\t{joined_ref_ids}')




if __name__ == "__main__":
    main()
    



