'''this module could parse the default bowtie output and organized into one .csv file, which use microRNA id as rowname, sample name as colname, and the panel in the table means the frequency of that kind of microRNA'''
'''Note: When encounter the region like 1)start point 100 & 106, 2)end point 120 & 126, and the overlaps from the sequence both satisfy conditions (100-120, 106-126),
the function Pendulumcheck will only count once instead of twice'''
from collections import defaultdict
import os
import pandas as pd

def PostionGet(posinput):
    '''this function is used to parse the position file which generated from miRBase hairpin and mature fasta (which acquired from .pos file)
    The position file contains 4 columns, and it is tab-delimitted.
    1) the ID of premature microRNA, which is indicated by the string 'MIR'
    2) the start postion of mature microRNA with respect to premature microRNA
    3) the end postion of mature microRNA with respect to premature microRNA
    4) The ID of mature microRNA, which is indicated by the string 'miR'. '''
    premature2start = defaultdict(list)
    premature2end = defaultdict(list)
    premature2mature = defaultdict(list)

    with open(posinput, 'rt') as input:
        for line in input:
            parsed_line = line.strip('\n').split('\t')
            tmp_key = parsed_line[0]
            tmp_start = parsed_line[1]
            tmp_end = parsed_line[2]
            tmp_mature = parsed_line[3]
            
            premature2start[tmp_key].append(tmp_start)
            premature2end[tmp_key].append(tmp_end)
            premature2mature[tmp_key].append(tmp_mature)
                
    input.close()
    return premature2start, premature2end, premature2mature

    
    # with open(posinput, 'rt') as input:
    #     for line in input:
    #         parsed_line = line.strip('\n').split('\t')
    #         premature2start[ parsed_line[0] ] = int(parsed_line[1])
    #         premature2end[ parsed_line[0] ] = int(parsed_line[2])
    #         premature2mature[ parsed_line[0] ] = parsed_line[3]
    # input.close()
    # return premature2start, premature2end, premature2mature

def Intervaljudge(l):
    '''the valid microRNA length is [18, 24], seq not in this interval will be assigned True flag'''
    lowerlimit = 18
    upperlimit = 24
    if l < 18 or l > 24:
        return True
    else:
        return False


def Pendulumcheck(alignment, premature2start, premature2end):
    '''this funciton have several conditions to check whether the alignent is valid or not.
    Note: the input of this function is parsed by Python split using Tab.
    1) end point of input microRNA minors the start point of standard sequence in the miRBase should be greater than 18
    2) the end point of standard sequence in the miRBase minors the start point of input microRNA should be greater than 18
    3) the input microRNA sequence is within the interval of standard miRBase sequence and its own lengths should be greater than 18 but smaller than 24'''
    # Save the input microRNA information
    microRNA_length = int(alignment[0].split('_')[2])
    # microRNA_start = int(alignment[3]) + 1
    microRNA_start = int(alignment[3])
    microRNA_end = microRNA_start + microRNA_length
    microRNA_count = int(alignment[0].split('_')[1])
    index = alignment[2]
    # print(microRNA_count)
    # print(f"len: {microRNA_length}, start: {microRNA_start}, end: {microRNA_end}, index: {index}")
    # Index the standard microRNA information
    '''Because Now I am using PimREN database, but the position file I am using is from miRBase
    there are some difference, like the name of the premature microRNA
    e.g. 
    PimREN: Osa-MIR168
    miRBase: osa-MIR168a, osa-MIR168b, which means there are more details in the miRBase
    '''
    SDmicroRNA_start = premature2start[index]
    SDmicroRNA_end = premature2end[index]
    # SDmicroRNA_start = int(premature2start[index])
    # SDmicroRNA_end = int(premature2end[index])
    # print(len(SDmicroRNA_start))
    # print(f"SDstart: {SDmicroRNA_start}, SDeNnd: {SDmicroRA_end}")

    # Check whether the hairpin has 5' + 3' or just one mature microRNA 
    if len(SDmicroRNA_start) == 1:
        '''Set a expression to find out whether it is a valid alignment or not
        version1'''
        # expression1 = ((microRNA_end - int(SDmicroRNA_start) >= 18) and (microRNA_end <= int(SDmicroRNA_end)))
        # expression2 = ((int(SDmicroRNA_end) - microRNA_start >= 18) and (microRNA_start >= int(SDmicroRNA_start)))
        # expression3 = (microRNA_start >= int(SDmicroRNA_start) and microRNA_end <= int(SDmicroRNA_end))
        # expression4 = (microRNA_start <= int(SDmicroRNA_start) and microRNA_end >= int(SDmicroRNA_end))
        '''Set a expression to find out whether it is a valid alignment or not
        version2'''
        expression1 = ((microRNA_end - int(SDmicroRNA_start[0]) >= 18))
        expression2 = ((int(SDmicroRNA_end[0]) - microRNA_start >= 18))
        if expression1 and expression2:
            return True, 0, microRNA_count
        else:
            return False, 0, microRNA_count
    
    elif len(SDmicroRNA_start) == 2:
        '''Set number for the 5' and 3' microRNA,
        if a hairpin has 5' and 3' mature microRNA, We give one the sequence with index number "1", another sequence with index number "2" '''
        SDmicroRNA_start_1 = int(SDmicroRNA_start[0])
        SDmicroRNA_start_2 = int(SDmicroRNA_start[1])
        SDmicroRNA_end_1 = int(SDmicroRNA_end[0])
        SDmicroRNA_end_2 = int(SDmicroRNA_end[1])

        '''Set a expression to find out whether it is a valid alignment or not
        version1'''
        # expression1 = ((microRNA_end - SDmicroRNA_start_1) >= 18 and (microRNA_end <= SDmicroRNA_end_1))
        # expression2 = ((SDmicroRNA_end_1 - microRNA_start) >= 18 and (microRNA_start >= SDmicroRNA_start_1))
        # expression3 = (microRNA_start >= SDmicroRNA_start_1 and microRNA_end <= SDmicroRNA_end_1)
        # expression4 = (microRNA_start <= SDmicroRNA_start_1 and microRNA_end >= SDmicroRNA_end_1)
        
        # expression5 = ((microRNA_end - SDmicroRNA_start_2) >= 18 and (microRNA_end <= SDmicroRNA_end_2))
        # expression6 = ((SDmicroRNA_end_2 - microRNA_start) >= 18 and (microRNA_start >= SDmicroRNA_start_2))
        # expression7 = (microRNA_start >= SDmicroRNA_start_2 and microRNA_end <= SDmicroRNA_end_2)
        # expression8 = (microRNA_start <= SDmicroRNA_start_2 and microRNA_end >= SDmicroRNA_end_2)

        '''Set a expression to find out whether it is a valid alignment or not
        version2'''
        expression1 = ((microRNA_end - SDmicroRNA_start_1) >= 18)
        expression2 = ((SDmicroRNA_end_1 - microRNA_start) >= 18)
        expression3 = ((microRNA_end - SDmicroRNA_start_2) >= 18)
        expression4 = ((SDmicroRNA_end_2 - microRNA_start) >= 18)
        if expression1 and expression2:
            return True, 0, microRNA_count
        elif expression3 and expression4:
            return True, 1, microRNA_count
        else:
            return False, 0, microRNA_count
    
    elif len(SDmicroRNA_start) == 3:
        SDmicroRNA_start_1 = int(SDmicroRNA_start[0])
        SDmicroRNA_start_2 = int(SDmicroRNA_start[1])
        SDmicroRNA_start_3 = int(SDmicroRNA_start[2])
        SDmicroRNA_end_1 = int(SDmicroRNA_end[0])
        SDmicroRNA_end_2 = int(SDmicroRNA_end[1])
        SDmicroRNA_end_3 = int(SDmicroRNA_end[2])

        expression1 = ((microRNA_end - SDmicroRNA_start_1) >= 18)
        expression2 = ((SDmicroRNA_end_1 - microRNA_start) >= 18)
        expression3 = ((microRNA_end - SDmicroRNA_start_2) >= 18)
        expression4 = ((SDmicroRNA_end_2 - microRNA_start) >= 18)
        expression5 = ((microRNA_end - SDmicroRNA_start_3) >= 18)
        expression6 = ((SDmicroRNA_end_3 - microRNA_start) >= 18)
        if expression1 and expression2:
            return True, 0, microRNA_count
        elif expression3 and expression4:
            return True, 1, microRNA_count
        elif expression5 and expression6:
            return True, 2, microRNA_count
        else:
            return False, 0, microRNA_count



def OrganizeBDoutput(bowtiedefault, premature2start, premature2end, premature2mature):
    '''this function is used to parse the bowtie default output so as to count the expression of each mature microRNA
    the main function of OrganizeBDoutput, is like below:
    1) filter the alignment which aligned reversely to the reference premature microRNA
    2) save the valid microRNA alignment'''
    four3b = open('444b.woziji.txt', 'w')
    maturemicroRNA_freq_dict = {}
    # unmatchmicroRNA_freq_dict = {}
    filter_count = 0
    unmatched_count = 0
    with open(bowtiedefault, 'rt') as input:
        # unvalid_count = 0
        for line in input:
            # print(int(line.split('\t')[0].split('_')[2]))
            # print(line.split('\t')[1])
            ''' 1. to check whether the alignment is forward or reverse,
                only save the forward alignment
                2. to filter the seq whose length is not in the interval of [18, 24]
                3. We allow 2 bp shift'''
            try:
                tmp_len = int(line.split('\t')[0].split('_')[2])
                tmp_icon = line.split('\t')[1]
                if  Intervaljudge(tmp_len)  or  tmp_icon == '-':
                    tmp_name = line.split('\t')[0]
                    print(f"{tmp_name} len: {tmp_len}, type: {type(tmp_len)}, icon: {tmp_icon}, this seq is too long or not forward aligned")
                    # unvalid_count += 1
                    continue
                flag, tmp_index, microRNA_count = Pendulumcheck(line.split('\t'), premature2start, premature2end)
                print(f"{flag}, {tmp_index}, {microRNA_count}")
                # print(f"{flag}, {tmp_index}")
                if flag:
                    '''this part need revision, because the some hairpin has two mature microRNA'''
                    index = line.split('\t')[2].split('-')[0].lower() + '-' + line.split('\t')[2].split('-')[1]
                    if index == "osa-MIR444b" and premature2mature[index][tmp_index] == "osa-miR444b.1":
                        four3b.write(line)
                    # print(f"{index}, {tmp_index}")
                    maturemicroRNA_freq_dict[premature2mature[index][tmp_index]] = maturemicroRNA_freq_dict.get(premature2mature[index][tmp_index], 0) + microRNA_count
                else:
                    filter_count += microRNA_count
            except:
                tmp_name = line.split('\t')[0]
                print(f"{tmp_name} could not be anchored")
                unmatched_count += 1
    print(f"[Warning] {unmatched_count} Unmatched")
    print(f"Filter {filter_count} sequences")
    return maturemicroRNA_freq_dict


def OneBowtieOutputOrganize(sname):
    tmp_dict = OrganizeBDoutput(sname, premature2start, premature2end, premature2mature)
    tmp_df = pd.DataFrame(data=tmp_dict.values(), index = list(tmp_dict))
    return tmp_df


def BowtieOutputNameGet(path):
    '''this function is used to grep all the filename under the dir bowtie_output'''
    tmp_list = []
    for x in os.listdir(path):
        tmp_list.append(path + '/' + x)
    return tmp_list


def MultipleBowtieOuputOrganize(inputlist):
    '''this function is used to run the OrganizeBDoutput multiple times (deal multiple samples at the same time and generate the output using pandas)'''
    tmp_df = pd.DataFrame()
    for x in inputlist:
        sample_name = x.split('/')[1].split('.')[0]
        tmp_dict = OrganizeBDoutput(x, premature2start, premature2end, premature2mature)
        tmp_df[sample_name] = pd.Series(tmp_dict.values(), index=tmp_dict.keys())
    
    return tmp_df



if __name__ == "__main__":
    
    premature2start, premature2end, premature2mature = PostionGet('osa_hairpin_mature.pos')

    bowtieoutput_list = BowtieOutputNameGet('bowtie_output')
    print(bowtieoutput_list)

    # tmp_df = OneBowtieOutputOrganize('bowtie_output/32RT-1.2207091646.sam')
    # print(tmp_df)
    # tmp_df.to_csv('./tmp.32RT-1.csv')

    tmp_df = MultipleBowtieOuputOrganize(bowtieoutput_list)
    tmp_df.to_csv('all_samples.csv')
    print(tmp_df)
