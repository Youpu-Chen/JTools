'''this module is used to parse the information from miRBase database
Here are the main goal of this module
1) parse the informtion in the organism text file
2) generathe the hairpin and mature RNA of specified species '''

import gzip

species_name = 'osa'
opener = gzip.open

def LatinNamematch(oginput):
    with opener(oginput, 'rt') as input:
        for line in input:
            if '#' in line:
                continue
            else:
                parsed_line = line.strip('\n').split('\t')
                if parsed_line[0] == species_name:
                    return parsed_line[2]
                    break
    input.close()


def Organismhairpin():
    '''this function generates the hairpin sequence of the species user specified
    Note: because the raw data we use is in DNA format, so we transform the sequence from '''
    seqid_list = []
    seq_list = []
    tmp_seq = ''
    with opener('hairpin.fa.gz', 'rt') as input:
        for index, line in enumerate(input):
            parsed_line = line.strip('\n')
            if ">" in parsed_line and index == 0:
                seqid_list.append(parsed_line)
            elif ">" in parsed_line and index != 0:
                seqid_list.append(parsed_line)
                seq_list.append(tmp_seq)
                tmp_seq = ''
            else:
                tmp_seq += parsed_line.replace('U', 'T')
        seq_list.append(tmp_seq)
    input.close()

    outputname = species_name + '.hairpin.fa'
    with open(outputname, 'w') as output:
        for index, id in enumerate(seqid_list):
            if latin_name in id:
                output.write(id + '\n')
                output.write(seq_list[index] + '\n')
    output.close()

def Organismmature():
    '''this function generates the mature microRNA sequence of the species user specified'''
    seqid_list = []
    seq_list = []
    tmp_seq = ''
    with opener('mature.fa.gz', 'rt') as input:
        for index, line in enumerate(input):
            parsed_line = line.strip('\n')
            if ">" in parsed_line and index == 0:
                seqid_list.append(parsed_line)
            elif ">" in parsed_line and index != 0:
                seqid_list.append(parsed_line)
                seq_list.append(tmp_seq)
                tmp_seq = ''
            else:
                tmp_seq += parsed_line.replace('U', 'T')
        seq_list.append(tmp_seq)
    input.close()

    outputname = species_name + '.mature.fa'
    with open(outputname, 'w') as output:
        for index, id in enumerate(seqid_list):
            if latin_name in id:
                output.write(id + '\n')
                output.write(seq_list[index] + '\n')
    output.close()


# def Getdead(deadinput):
#     '''this function is useless, because the data download from miRBase has already removed them'''
#     seqid_list = []
#     # seq_list = []
#     with opener(deadinput, 'rt') as input:
#         for line in input:
#             tmp_line = line.strip('\n')
#             # print(tmp_line)
#             if 'ID' in tmp_line and species_name in tmp_line:
#                 seqid_list.append(tmp_line.split(' ')[3])
#     input.close()
#     return seqid_list




if __name__ == "__main__":

    latin_name = LatinNamematch('organisms.txt.gz')
    Organismhairpin()
    Organismmature()
