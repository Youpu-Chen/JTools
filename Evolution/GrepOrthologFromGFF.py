'''this module use rewrite Ortholog file to grep corresponding ID'''
import sys


def ParseOrthologs(filename):
    geneids = []
    with open(filename, 'rt') as input:
        for line in input:
            geneids.append(line.strip('\n'))
    return geneids


def ParseGeneGFF3(geneids, filename, outputfile):
    output = open(outputfile, 'w')
    with open(filename, 'rt') as input:
        for line in input:
            parsed_line = line.strip('\n').split('\t')
            id = parsed_line[len(parsed_line) - 1]
            # print(id)
            if id.strip('ID=') in geneids:
                print(id)
                chr = parsed_line[0]
                start = parsed_line[3]
                end = parsed_line[4]
                info = parsed_line[5]
                strand = parsed_line[6]
                output.write(f'{chr}\t{start}\t{end}\t{id.strip("ID=")}\t{info}\t{strand}\n')
    output.close()


geneids = ParseOrthologs(sys.argv[1])  # specify rewrite.Orthogroups.tsv
print(geneids[:5])
ParseGeneGFF3(geneids, sys.argv[2], sys.argv[3])  # specify gff3 file[2] and output file name [3]




                