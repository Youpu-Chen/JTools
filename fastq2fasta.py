#---------------------------------------------------------------------------------------------
# Description
# this script is used to transfer single-ended fastq to fasta, and 
# reformat the sequence id using frequency and 
#---------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
# Data Preparation
# 1) *.fq or *.fq.gz
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
# Test Code
#---------------------------------------------------------------------------------------------
# uniq_seq = []
# seq_count = {}
# with open('test-32RT-1.fq.gz', 'r') as input:
    
#     for line_index, line in enumerate(input):

#         if line_index % 4 == 1:
#             tmp_line = line.strip('\n')
#             # seq_list.append(tmp_line)
#             seq_count[tmp_line]= seq_count.get(tmp_line, 0) + 1
            
#         else:
#             continue

# uniq_seq = seq_count.keys()
# input.close()




#---------
# set args
#---------
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gzip", action="store_true", required=False, help="Set if the input fastq is gzipped.")
parser.add_argument('-i', "--input", type=str, metavar='Input fastq file', required=True, help='It should be the single-ended Illumina sequence')
parser.add_argument('-o', "--output", type=str, metavar='Output fasta', required=True, help='Output name of fasta file')
args = parser.parse_args()

if args.gzip:
    opener = gzip.open
else:
    opener = open


#---------
# analysis
#---------
def CountandFormattingseq(file_name, output_name):
    ### formatting
    uniq_seq = []
    seq_count = {}
    with opener(file_name, 'rt') as input:
        for line_index, line in enumerate(input):
            if line_index % 4 == 1:
                tmp_line = line.strip('\n')
                seq_count[tmp_line]= seq_count.get(tmp_line, 0) + 1
            else:
                continue
    uniq_seq = seq_count.keys()
    input.close()
    
    ### write 
    with open(output_name, 'w') as output:
        for index, seq in enumerate(uniq_seq, 1):
            output.write(f'>seq_{index}_{seq_count[seq]}_{len(seq)}\n')
            output.write(f'{seq}\n')
    output.close()

    return print("The script run successfully!")

CountandFormattingseq(args.input, args.output)
