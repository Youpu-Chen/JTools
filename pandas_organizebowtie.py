import enum
from turtle import Turtle
from unicodedata import numeric
from organizebowtie import BowtieOutputNameGet
import pandas as pd

def ParsePos(pname):
    tmp_input = pd.read_table(pname, sep='\t', names=['hpID', 'micro_start', 'micro_end', 'microID'])
    df1 = pd.DataFrame(tmp_input)
    # print(df1)
    return df1



def OrganizeBowtieOutput(bname, p_df):
    tmp_input = pd.read_table(bname, sep='\t', names=['seqid', 'direction', 'hpID', 'seq_start', 'seq', 'A', 'B', 'C', 'D'])
    df1 = pd.DataFrame(tmp_input)
    # Only select columns with usable columns
    df2 = df1.loc[:, ['seqid', 'direction', 'hpID', 'seq_start']]

    # Merge the sample data using the key hpID
    df3 = pd.merge(df2, p_df, on='hpID')
    # print(df3)
    # df3.to_csv('df3_test.csv')
    df3['seq_length'] = pd.to_numeric(df3['seqid'].str.split('_', expand=True)[2])
    df3['count'] = pd.to_numeric(df3['seqid'].str.split('_', expand=True)[1])
    # df3.to_csv('df3_test.csv')
    df3['seq_end'] = df3['seq_start'] + df3['seq_length']
    # df3.to_csv('df3_test.csv')

    # Using the conditions to filter the unvalid alignment 
    df3['overlap1'] = df3['seq_end'] - df3['micro_start']
    df3['overlap2'] = df3['micro_end'] - df3['seq_start']
    # df3.to_csv('df3_test.csv')

    # Filter raw data using condition above
    df4 = df3[(df3['overlap1'] >= 18) & (df3['overlap2'] >= 18) & (df3['direction'] != "-")]
    df5 = df4[((df4['seq_length'] >= 18) & (df4['seq_length'] <= 24))]
    # print(df5)
    # df5.to_csv('df5_test.csv')
    
    # # Group the filter data and do the summary statistics
    # micro_count = df5.groupby('microID')['count'].sum()  # if not using pd.DataFrame(), this function will return a Series
    micro_count = pd.DataFrame(df5.groupby('microID')['count'].sum())
    # print(micro_count.index.name)
    # print(micro_count.name)
    # micro_count.to_csv('micro_test.csv')
    return micro_count


def MultipleBowtieOuputOrganize(inputlist, pos):
    '''the inputlist is the list which contains the name of bowtie output'''
    # 1. Using pd.concat()
    tmp_all_df = pd.DataFrame()
    for x in inputlist:
        tmp_df = OrganizeBowtieOutput(x, pos)
        tmp_df.columns = [x.split('/')[1].split('.')[0]]
        tmp_all_df = pd.concat([tmp_all_df, tmp_df], axis=1)
    tmp_all_df.to_csv("pandas_all_samples.csv")

    # 2. Using merge()
    # tmp_all_df = pd.DataFrame()
    # # print(tmp_all_df)
    # for index, x in enumerate(inputlist):
    #     if index == 0:
    #         tmp_df = OrganizeBowtieOutput(x, pos)
    #         # print(tmp_df.index.name)
    #         tmp_all_df[x.split('/')[1].split('.')[0]] = tmp_df
    #         # print(tmp_all_df)
    #         # print(tmp_all_df.index.name)
    #     else:
    #         tmp_df = OrganizeBowtieOutput(x, pos)
    #         tmp_df.columns = [x.split('/')[1].split('.')[0]]
    #         # print(tmp_df.columns)
    #         # print(tmp_df)
    #         tmp_all_df = pd.merge(tmp_all_df, tmp_df, left_index=True, right_index=True, how='outer')
    #         # print(tmp_all_df)   
    # # print(tmp_all_df)
    # tmp_all_df.to_csv("pandas_all_samples.csv")

    # 3. Using join()
    # tmp_all_df = pd.DataFrame()
    # for index, x in enumerate(inputlist):
    #     if index == 0:
    #         tmp_df = OrganizeBowtieOutput(x, pos)
    #         tmp_all_df[x.split('/')[1].split('.')[0]] = tmp_df
    #     else:
    #         tmp_df = OrganizeBowtieOutput(x, pos)
    #         tmp_df.columns = [x.split('/')[1].split('.')[0]]
    #         # print(tmp_df.columns)
    #         tmp_all_df = tmp_all_df.join(tmp_df, how='outer')
    # tmp_all_df.to_csv("pandas_all_samples.csv")

    # 4. Using merge() and not using loop
    '''this version will report future warning'''
    # tmp_all_df = pd.DataFrame()
    # # print(tmp_all_df)
    # for x in inputlist:
    #     tmp_df = OrganizeBowtieOutput(x, pos)
    #     tmp_all_df = pd.merge(tmp_all_df, tmp_df, left_index=True, right_index=True, how='outer')
    #         # print(tmp_all_df)   
    # # print(tmp_all_df)
    # tmp_all_df.to_csv("pandas_all_samples.csv")






if __name__ == "__main__":
    pos_df = ParsePos('osa_hairpin_mature.pos')
    # OrganizeBowtieOutput('bowtie_output/32RT-1.2207091646.sam', pos_df)
    bowtieoutput_list = BowtieOutputNameGet('bowtie_output')
    MultipleBowtieOuputOrganize(bowtieoutput_list, pos=pos_df)