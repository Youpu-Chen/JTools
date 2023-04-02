import pandas as pd
pd.options.display.max_columns = None

import seaborn as sns
import matplotlib.pyplot as plt


def readFreq(input_file):
    freq_list = []
    input = open(input_file, "rt")
    for line in input:
        parsed_line = line.strip("\n").split("\t")
        # print(line)

        # 过滤sex chromosome上的位点
        if parsed_line[0] == "chrX" or parsed_line[0] == "chrY":
            continue

        # 只保留biallelic SNP
        if len(parsed_line[ len(parsed_line) - 1 ].split(";")) > 2:
            continue
        
        # 过滤INDEL
        for allele_records in parsed_line[ len(parsed_line) - 1 ].split(";"):
            allele_type = allele_records.split(":")[0]
            if len(allele_type) > 1:
                break

        # 拆分得到REF和ALT allele
        # parsed_line[-1] = parsed_line[-1].split(";")
        for allele_type in parsed_line[-1].split(";"):
            parsed_line.append(allele_type)
        # print(parsed_line)

        # 对通过筛选条件的行添加tag信息
        tag = parsed_line[0] + "-" + str(parsed_line[1])
        parsed_line.append(tag)
        # print(parsed_line)

        # 去除原始的allele freq信息
        del parsed_line[-4]

        # print(parsed_line)
        freq_list.append(parsed_line)

    return freq_list

EUR_list = readFreq("./concat_freq_output/concat.EUR.all.AF.frq")
EAS_list = readFreq("./concat_freq_output/concat.EAS.all.AF.frq")


def readFstSite(input_file):
    fst_df = pd.read_table(input_file, sep="\t", header=None)
    fst_df.columns = ["Chr", "Start", "End", "N_VARIANTS", "Weighted_Fst", "Mean_Fst"]
    # print(fst_df.head())

    # 索引top 5% Fst
    top5_fst_df = fst_df[ fst_df["Weighted_Fst"] >= fst_df["Weighted_Fst"].quantile(0.95) ]
    top5_fst_df["Fst_Tag"] = "HD-Fst"

    # 添加Pos信息
    top5_fst_df["Tag"] = top5_fst_df["Chr"] + "-" + top5_fst_df["Start"].astype(str)

    # 重新整合得到top5_fst_df
    reorder_top5_fst_df = top5_fst_df[["Tag", "Fst_Tag", "Fst_Tag", "Weighted_Fst"]]

    # print(reorder_top5_fst_df.head())
    return reorder_top5_fst_df

reorder_top5_fst_df = readFstSite("./EUR.vs.EAS.all.fst")


def CompareAndTag(pop1_list, pop2_list, fst_df):
    pop1_df = pd.DataFrame(pop1_list)
    pop2_df = pd.DataFrame(pop2_list)

    # 给pop df命名
    pop1_df.columns = ["Chr", "Start", "N", "S", "REF_Allele", "ALT_Allele", "Tag"]
    pop2_df.columns = ["Chr", "Start", "N", "S", "REF_Allele", "ALT_Allele", "Tag"]

    # 
    freq_df = pop1_df[["Tag", "REF_Allele", "ALT_Allele"]].merge(pop2_df[["Tag", "REF_Allele", "ALT_Allele"]], on="Tag", how="outer")


    # 只保留HD-Fst sites
    HD_Fst_df = fst_df.merge(freq_df, on="Tag", how="left")
    # print(HD_Fst_df.head())

    # 拆分得到allele frequency
    HD_Fst_df["REF_Allele_x"] = HD_Fst_df["REF_Allele_x"].str.split(":").str[1].astype(float)
    HD_Fst_df["ALT_Allele_x"] = HD_Fst_df["ALT_Allele_x"].str.split(":").str[1].astype(float)
    HD_Fst_df["REF_Allele_y"] = HD_Fst_df["REF_Allele_y"].str.split(":").str[1].astype(float)
    HD_Fst_df["ALT_Allele_y"] = HD_Fst_df["ALT_Allele_y"].str.split(":").str[1].astype(float)
    
    print(HD_Fst_df.head())
    return HD_Fst_df


HD_Fst_df = CompareAndTag(EUR_list, EAS_list, reorder_top5_fst_df)
HD_Fst_df.to_csv("./HD_Fst_df", sep="\t", index=False)

def createDensityPlot(input_df):
    # set a grey background (use sns.set_theme() if seaborn version 0.11.0 or above) 
    sns.set(style="darkgrid")

    # plotting both distibutions on the same figure
    fig = sns.kdeplot(input_df['REF_Allele_x'], shade=True, color="red")
    fig = sns.kdeplot(input_df['REF_Allele_y'], shade=True, color="blue")
    fig = sns.kdeplot(input_df['ALT_Allele_x'], shade=True, color="yellow")
    fig = sns.kdeplot(input_df['ALT_Allele_y'], shade=True, color="green")
    density_plot = fig.get_figure()
    density_plot.savefig("./EUR-EAS-allele-frequency-highly-differentiated.pdf", format="pdf", bbox_inches="tight")
    # plt.show()


createDensityPlot(HD_Fst_df)