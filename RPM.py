'''this module is used to convert the raw count table to RPM normalized table'''
import pandas as pd
import os



def Countsum(bwtname):
    # with open(bwtname, 'rt') as input:
    bwt_df = pd.DataFrame(pd.read_table(bwtname, index_col=False, names=['seqID', 'direction', 'refID', 'start_pos', 'seq', 'quality', 'randomN']))
    return pd.to_numeric(bwt_df['seqID'].str.split('_', expand=True)[1]).sum()


def BowtieOutputNameGet(path):
    '''Note: This function is from organizebowtie.py
    this function is used to grep all the filename under the dir bowtie_output'''
    tmp_list = []
    for x in os.listdir(path):
        tmp_list.append(path + '/' + x)
    return tmp_list


def ColRPMconvert(x, count):
    '''this function is applied to one column of pandas table
    x means the column variable in DataFrame
    count is the value returned by the function Countsum'''
    return (x * 1000000 / count)


def DFRPMconvert(rawcount, bwtlist):
    # Read the raw count table
    rawcount_df = pd.DataFrame(pd.read_csv(rawcount, index_col=0, header=0))
    rawcount_df.fillna(0, inplace=True)
    # print(rawcount_df)
    # print(rawcount_df['32RT-1'].apply(ColRPMconvert, args=(count,)))

    # Using loop to apply the ColRPMconvert to each raw column
    for x in bwtlist:
        bwt_name = x.split('/')[1].split('.')[0]
        count = Countsum(x)
        rawcount_df[bwt_name] = rawcount_df[bwt_name].apply(ColRPMconvert, args=(count,))
    
    rawcount_df.to_csv('RPM.normalized.csv')



if __name__ == "__main__":
    # tmp_sum = Countsum('test.sam')
    # print(type(tmp_sum))
    # tmp_sum = Countsum('32RT-1.2206301952.sam')
    # print(tmp_sum)
    
    bwtlist = BowtieOutputNameGet('bowtie_output')
    # print(bwtlist)

    DFRPMconvert('pandas_all_samples.csv', bwtlist)