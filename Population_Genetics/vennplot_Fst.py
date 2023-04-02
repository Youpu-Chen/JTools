'''
Description: Venn plot for pair-wise comparison of Fst
Note:
- manually specify pyvenn path
'''
pyvenn_path="./pyvenn"

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import rcParams
import sys
sys.path.append(pyvenn_path)
import venn

print("pakages load successfully!")

# basis set
rcParams['font.family'] = 'Microsoft YaHei'

def readFst(filename):
    '''read Fst with each pair of Fst calculation concatenated from result of vcftools'''
    Fst = pd.read_table(filename, sep="\t", 
                        names=["Chr", "Start", "End", "N_VARIANTS", "Weighted_Fst", "Mean_Fst"],
                        dtype={"Chr": str,
                               "Start": int,
                               "End": int,
                               "N_VARIANTS": int,
                               "Weighted_Fst": np.float64, 
                               "Mean_Fst": np.float64
                               },
                        header=None)  # ignore header

    # print(Fst.dtypes)
    # print(Fst.head(5))
    return Fst

def filterFst(dataframe):
    '''only saved top 5% of pair-wise Fst values'''
    dataframe = dataframe[ dataframe["Weighted_Fst"] >= dataframe["Weighted_Fst"].quantile(0.95) ]
    return dataframe

def addPos(dataframe):
    '''apply "chromosome-pos" as the tag for further analysi
    '''
    # 添加Pos信息
    dataframe["Pos"] = dataframe["Chr"] + "-" + dataframe["Start"].astype(str)  # make sure the datatype is "str"
    # print(dataframe.head())
    return dataframe


def vennPlot(dataframes, output_picture, tag):
    '''
    Note: dataframe is a list which contains multiple dataframe for Venn plot
    '''
    labels = venn.get_labels(dataframes, 
                             fill=["number", "percent"])
    fig, ax = venn.venn5(labels, names=tag)
    fig.savefig(output_picture, format="pdf", bbox_inches="tight")
    # fig.show()

def parseTag(filename):
    '''parse the input name of specified fiel
    e.g., FIN.vs.KHV.all.fst, the target is KHV.
    '''
    print(filename.split(".")[2])
    return filename.split(".")[2]

def main():
    '''parse input from command lines'''
    # load dataset
    fn_list = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]]  # default for 5-D venn plot
    Fst_list = []
    Tag_list = []
    for file in fn_list:
        Fst = readFst(file)
        Tag = parseTag(file)
        
        # Fst filtering
        filtered_Fst = filterFst(Fst)
        tag_filter_Fst = addPos(filtered_Fst)
        
        Tag_list.append(Tag)
        Fst_list.append(set(tag_filter_Fst["Pos"]))
    print("load dataset successfully!")

    # make plot
    vennPlot(Fst_list, sys.argv[6], Tag_list)
    print("plot venn, successfully!")
    
    

if __name__ == "__main__":
    main()