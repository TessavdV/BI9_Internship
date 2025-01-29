import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

"""
CADDSV FT and CT comparison
"""

def amount_scored(input_files, samplename, output, score):
    """
    Creates barplot showing total variants, amount scored by FT and by CT.
    Not-scored seperated in located on Y-chr and not located on Y-chr.
    """
    try:
        combined_data = []

        # Loop through input files to read and combine data
        for file in input_files:
            df = pd.read_csv(file, sep='\t')
            combined_data.append(df)

        # Concatenate all dataframes
        combined_df = pd.concat(combined_data, ignore_index=True)

        # Extract relevant columns
        ft_scores = combined_df.iloc[:, 8]  # Column 9 ('CADDSV_FT_Score')
        if score == 'CT-O':
            score_column = 11 # Column 12 ('MAX_OVERLAP')

        elif score == 'CT-P':
            score_column = 9 # Column 10 ('MAX_PATH_SCORE')

        ct_scores = combined_df.iloc[:, score_column]

        # Count total, FT-scored, and CT-scored variants
        total_amount = len(combined_df)
        ft_count = ft_scores.notnull().sum()
        ct_count = ct_scores.notnull().sum()

        # Count Y-chromosome null scores
        subset_ychr = combined_df[combined_df["CHROM"] == "Y"]
        y_chr_ft_null = subset_ychr.iloc[:, 8].isnull().sum()
        y_chr_ct_null = subset_ychr.iloc[:, score_column].isnull().sum()

        # Count nonY-chromosome null scores
        subset_ychr = combined_df[combined_df["CHROM"] != "Y"]
        noy_chr_ft_null = subset_ychr.iloc[:, 8].isnull().sum()
        noy_chr_ct_null = subset_ychr.iloc[:, score_column].isnull().sum()

        # Data for the bar plot
        labels = ['Total Variants', 'FT-Scored', 'CT-Scored']
        bar_values = np.array([total_amount, ft_count, ct_count])
        stacked_values = np.array([0, y_chr_ft_null, y_chr_ct_null])
        stacked_twovalues = np.array([0, noy_chr_ft_null, noy_chr_ct_null])

        # Create the bar plot
        plt.figure(figsize=(8, 6))
        plt.bar(labels, bar_values, color='steelblue', label='Scored')
        plt.bar(labels, stacked_values, bottom=bar_values, color='darkorange', label='Not scored from Y-chr')
        plt.bar(labels, stacked_twovalues, bottom=bar_values + stacked_values, color='mediumpurple',
                label='Not scored non-Y chr')

        # Add title, labels, and legend
        plt.title(samplename + ' ' + score + ': Variants Scored', fontsize=16)
        plt.ylabel('Number of Variants')
        plt.legend(loc='upper left')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()

        # Save and show the plot
        plt.savefig(output)
        plt.show()

        # Print amount of variants scored
        print(ct_count)

    except Exception as e:
        print(f"An error occurred: {e}")



if __name__ == "__main__":
    samplename_hifi = 'HifiCNV P4-C4'
    samplename_pbsv = "PBSV P13-H4"
    score_o = 'CT-O'
    score_p = 'CT-P'

    input_files_HifiCNV = [
        "C:/Users/Tessa.vanderVeer@radboudumc.nl/Documents/Result_analysis/CADDSV_FTvsCT/P4-C4_CADDSV_CTCPandFT.tsv",
    ]
    input_files_PBSV = [
        "C:/Users/Tessa.vanderVeer@radboudumc.nl/Documents/Result_analysis/CADDSV_FTvsCT/P13-H4_CADDSV_CTCPandFT.tsv",
    ]

    output_hifi_cto = "C:/Users/Tessa.vanderVeer@radboudumc.nl/Documents/Result_analysis/CADDSV_FTvsCT/results/"+samplename_hifi+'_'+score_o+"_CADDSV_FTvsCT.png"
    output_hifi_ctp = "C:/Users/Tessa.vanderVeer@radboudumc.nl/Documents/Result_analysis/CADDSV_FTvsCT/results/"+samplename_hifi+'_'+score_p+"_CADDSV_FTvsCT.png"
    output_pbsv_cto = "C:/Users/Tessa.vanderVeer@radboudumc.nl/Documents/Result_analysis/CADDSV_FTvsCT/results/"+samplename_pbsv+'_'+score_o+"_CADDSV_FTvsCT.png"
    output_pbsv_ctp = "C:/Users/Tessa.vanderVeer@radboudumc.nl/Documents/Result_analysis/CADDSV_FTvsCT/results/"+samplename_pbsv+'_'+score_p+"_CADDSV_FTvsCT.png"

    amount_scored(input_files_HifiCNV, samplename_hifi, output_hifi_cto, score_o)
    amount_scored(input_files_HifiCNV, samplename_hifi, output_hifi_ctp, score_p)
    amount_scored(input_files_PBSV, samplename_pbsv, output_pbsv_cto, score_o)
    amount_scored(input_files_PBSV, samplename_pbsv, output_pbsv_ctp, score_p)


