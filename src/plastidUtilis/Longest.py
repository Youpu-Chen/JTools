'''
Description: filter the unvaluable tRNA sequences of given .fa file
Note: "-d <delimiter>" parameter is used to retrieve the targeted ID instead of the whole header
'''
import optparse
usage = """python -m plastidUtilis.Longest -f <input.fasta> -d <delimiter> -o <output.fasta>
                                                                                    --Joe"""

parser = optparse.OptionParser(usage)
parser.add_option("-f", dest="input_fasta", help="input fasta",
                  metavar="FILE", action="store", type="string")
parser.add_option("-d", dest="delimiter", help="delimiter",
                  metavar="OUT", action="store", type="string")
parser.add_option("-o", dest="output", help="output fasta",
                  metavar="OUT", action="store", type="string")
(options, args) = parser.parse_args()

file_name = options.input_fasta
delimiter = str(options.delimiter)
output = options.output

fa_dict = {}
with open(file_name) as f:
    for line in f.readlines():
        if line.startswith(">"):
            seqname = line.strip()[1:]
            fa_dict[seqname] = []
        else:
            fa_dict[seqname].append(line.strip())
for k, v in fa_dict.items():
    fa_dict[k] = "".join(v)

keep_list = {}
for k in fa_dict.keys():
    gene = k.split(delimiter)[0]
    length = int(len(fa_dict[k]))
    if gene not in keep_list.keys():
        keep_list[gene] = k
    elif gene in keep_list and length > len(fa_dict[keep_list[gene]]):
        keep_list[gene] = k

with open(output, "w") as out:
    for v in keep_list.values():
        out.write(">" + v + "\n")
        out.write(fa_dict[v] + "\n")
out.close()