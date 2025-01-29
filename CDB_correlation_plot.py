"""
Creates CDB-test boxplot (incl. sex-chr data)
"""

import pandas as pd
import matplotlib.pyplot as plt


def create_boxplot(input_file, output_file, score):
    """
    Reads input file and creates boxplot.

    :param input_file: processed results of all CDB variants: AnnotSV, CT-O and CT-P score (and general info)
    :param output_file: boxplot file
    :param score: for which score boxplot is made (AnnotSV, CT-O or CT-P, or ACMG class)
    """
    df = pd.read_csv(input_file, sep='\t', header=0)

    required_columns = [score, 'CDB_CLASS']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Input file must contain the following columns: {', '.join(required_columns)}")

    # Haal INSERTION uit dataframe als SCORE = ANNOTSV
    if score == 'ANNOTSV':
        df = df[df['TYPE'] != 'INSERTION']

    # print hoeveel varianten geen score hebben gekregen
    NaDf_count = df[score].isna().sum()
    print("Amount of CDB variants without "+score+":  "+str(NaDf_count))

    df = df.dropna(subset=[score])

    class1_df = df.loc[df['CDB_CLASS'] == 'class 1']
    class2_df = df.loc[df['CDB_CLASS'] == 'class 2']
    class3_df = df.loc[df['CDB_CLASS'] == 'class 3']
    class4_df = df.loc[df['CDB_CLASS'] == 'class 4']
    class5_df = df.loc[df['CDB_CLASS'] == 'class 5']

    data1 = class1_df[score].tolist()
    data2 = class2_df[score].tolist()
    data3 = class3_df[score].tolist()
    data4 = class4_df[score].tolist()
    data5 = class5_df[score].tolist()
    fulldata = [data1, data2, data3, data4, data5]

    df_printje = class3_df[['CHROM','START', 'END', 'TYPE', 'ANNOTSV_SCORE', 'ANNOTSV_ACMG_CLASS', 'CDB_CLASS', 'MAX_PATH_SCORE', 'MAX_OVERLAP_SCORE']]
    df_printje = df_printje.sort_values(score)
    print(df_printje.to_string(index=False))

    plt.rcParams.update({'font.size': 14})

    # arrays doorgeven aan boxplot
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111)
    plt.boxplot(fulldata,
                medianprops=dict(lw=3, color='hotpink'))

    if score == 'ANNOTSV_SCORE' :
        plt.ylim(-2, 2)
        custom_ticks = [-1.5, -0.99, -0.9, 0, 0.9, 0.99, 1.5]
        plt.yticks(custom_ticks)
        plt.axhline(0.99, color='r', lw=1, linestyle='--')
        plt.axhline(-0.99, color='r', lw=1, linestyle='--')
        plt.axhline(0.90, color='C1', lw=1, linestyle=':')
        plt.axhline(-0.90, color='C1', lw=1, linestyle=':')
    elif score == 'ANNOTSV_ACMG_CLASS':
        plt.ylim(0.5, 5.5)
    else:
        plt.ylim(-0.5, 50)

    ax.set_xticklabels(['Class 1\n'+str(len(data1))+' vars', 'Class 2\n'+str(len(data2))+' vars',
                        'Class 3\n'+str(len(data3))+' vars', 'Class 4\n'+str(len(data4))+' vars',
                        'Class 5\n'+str(len(data5))+' vars'])
    plt.ylabel(score +' scores')
    plt.xlabel('CDB classes')
    plt.title('Correlation '+score+' scores and CDB-classification')


    plt.savefig(output_file)
    plt.show()

if __name__=="__main__":

    input_file = "C:/Users/Tessa.vanderVeer@radboudumc.nl/Documents/Result_analysis/CDB/CDB_fullresults_processed.tsv"
    scores = ['ANNOTSV_SCORE', 'MAX_PATH_SCORE', 'MAX_OVERLAP_SCORE', 'ANNOTSV_ACMG_CLASS']

    for score in scores:
        output_file="C:/Users/Tessa.vanderVeer@radboudumc.nl/Documents/Result_analysis/CDB/CDB_" + score + "_output.png"
        create_boxplot(input_file, output_file, score)