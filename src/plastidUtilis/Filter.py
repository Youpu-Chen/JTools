'''
Description: filter the unvaluable tRNA sequences of given .fa file
Note: set "-i 64" in general case.
'''
import optparse

usage = """python -m plastidUtilis.Filter -f <input.fasta> -i <minumum_length> -I <maximum_length> -o <output filename>
                                                                                                                    --Joe"""

parser = optparse.OptionParser(usage)
parser.add_option("-f", dest="input_file", help="fasta",
                  metavar="FILE", action="store", type="string")
parser.add_option("-i", dest="min", help="minimum length",
                  metavar="STRING", action="store", type="string")
parser.add_option("-I", dest="max", help="maximum length",
                  metavar="STRING", action="store", type="string")
parser.add_option("-o", dest="output", help="output file",
                  metavar="FILE", action="store", type="string")
(options, args) = parser.parse_args()

fa_dict = {}
fa_dict["dump"] = []
with open(options.input_file) as f:
    for line in f.readlines():
        if line.startswith(">"):
            header = line.strip()
            if header not in fa_dict.keys():
                fa_dict[header] = []
            elif header in fa_dict.keys():
                header = "dump"
        else:
            fa_dict[header].append(line.strip())
del fa_dict["dump"] 

with open(options.output, "w") as out:
    for k, v in fa_dict.items():
        fa_dict[k] = "".join(v)
        if len(fa_dict[k]) > int(options.min) and len(fa_dict[k]) < int(options.max):
            out.write(k + "\n")
            out.write(fa_dict[k] + "\n")
out.close()
                
        
        
