#---------------------------------------------------------------------------------------------
# Description
# This script is used for the addition of rgb color in JCVI result .simple file
#---------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
# Data Preparation
# 1) *.uniq.bed from JCVI
# 2) *.simple from JCVI
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
# Code
#---------------------------------------------------------------------------------------------
# gene_chr_list = []
# with open('moso.uniq.bed', 'r') as bed:
    
#     for line in bed.readlines():
#         tmp_line = line.strip('\n').split('\t')
#         gene = tmp_line[3]
#         chr = tmp_line[0]
#         link = gene + "\t" + chr + "\n"
#         gene_chr_list.append(link)

# with open('gene2chr.txt', 'w') as output:
#     output.writelines(gene_chr_list)

#---------
# set args
#---------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-b', "--bed", type=str, metavar='BED', required=True, help='BED file generated by JCVI')
parser.add_argument('-c', "--config", type=str, metavar='Color config', required=True, help='Personal defined RGB color configuration')
parser.add_argument('-s', "--simple", type=str, metavar='.simple', required=True, help='.simple file generated by JCVI')
parser.add_argument('-o', "--output", type=str, metavar='new.simple', required=True, help='Output name of new simple file')
args = parser.parse_args()


#---------
# acquire the gene position information
#---------
gene_chr_dic = {}
with open(args.bed, 'r') as bed:
    
    for line in bed.readlines():
        tmp_line = line.strip('\n').split('\t')
        gene = tmp_line[3]
        chr = tmp_line[0]
        gene_chr_dic[gene] = chr

bed.close


#---------
# set color config
#---------
color_dict = {}
with open(args.config, 'r') as color_input:
    for line in color_input.readlines():
        tmp_line = line.strip('\n').split('\t')
        pos = tmp_line[0]
        color_code = tmp_line[1]
        color_dict[pos] = color_code
color_input.close()



#---------
# add the color parameters to .simple file
#---------

with open(args.simple, 'r') as simple:
    with open(args.output, 'w') as output:
        for line in simple.readlines():
            stop_gene = line.split('\t')[1]
            pos = gene_chr_dic[stop_gene]
            color_set = color_dict[pos]

            output.writelines(color_set+"*"+line)
    output.close()
simple.close()