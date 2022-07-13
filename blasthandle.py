'''this module is used to handle the blast output, and generate the .pos file which is used to organize the results of bowtie alignment
Note: .pos file, the expecting format is like:
1st, the name of the premature miRNA sequence
2nd, the start position of the mature miRNA
3rd, the end position of the mature miRNA
4th, the ID of of mature miRNA sequence'''
import re


def Blastfilter(blastinput):
    '''This function is used to filter the blast default output "m6" based on the conditions below
    ) identity 100%
    ) the length of query sequence should within the interval of 20 and 24
    ) the number of gap and mismatches should be 0
    
    
    Each columns in fmt6 means: 1)qseqid 2)sseqid 3)pident 4)length 5)mismatch 6)gapopen 7)qstart 8)qend 9)sstart 10)send 11)evalue 12)bitscore
    '''

    with open('Filtered_' + blastinput, 'w') as output:
        with open(blastinput) as input:
            for line in input:
                parsed_line = line.strip('\n').split('\t')
                pident = parsed_line[2]
                milen = int(parsed_line[3])
                mismatch = parsed_line[4]
                gapopen = parsed_line[5]
                if mismatch != 0 and gapopen != 0 and pident != 100.000 and milen >= 20 and milen <= 24:
                    output.write(line)
    input.close()
    output.close()


def mirnaMatch(Fblastoutput):
    '''this function used the regular expression to match the microRNA info'''
    pattern1 = re.compile(r'([0-9]+[a-z]?)')
    pattern2 = re.compile(r'([0-9]+[a-z]?)')
    
    with open('Matched_' + 'Filtered_osa_mature_hairpin_blast.txt', 'w') as output:
        with open(Fblastoutput, 'r') as input:
            for line in input:
                parsed_line = line.rstrip().split("\t")
                mature = parsed_line[0]
                hairpin = parsed_line[1]
                match1 = pattern1.search(mature).group()
                match2 = pattern2.search(hairpin).group()
                if match1 == match2:
                    # print(line.rstrip())
                    output.write(line)


def Getpos(sname, MFblastoutput):
    '''this function generate the .pos file of specifed species, e.g. osa/rice'''
    output = open(sname + '_' + 'hairpin_mature.pos', 'w')
    with open(MFblastoutput, 'r') as input:
        for line in input:
            parsed_line = line.strip('\n').split('\t')
            hairpin = parsed_line[1]
            mature = parsed_line[0]
            startpos = parsed_line[8]
            endpos = parsed_line[9]
            output_line = '\t'.join([hairpin, startpos, endpos, mature]) + '\n'
            output.write(output_line)
    input.close()
    output.close()



if __name__ == "__main__":
    Blastfilter('osa_mature_hairpin_blast.txt')
    mirnaMatch('Filtered_osa_mature_hairpin_blast.txt')
    Getpos('osa', 'Matched_Filtered_osa_mature_hairpin_blast.txt')

