import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

"""
Plots accuracy algorithms for all files in dir. Calculates accuracy percentages. 
"""

def read_and_combine_files(input_dir, score):
    """
    Reads all the files in the given directory and combines them into one large DataFrame.
    Only files with the correct format (having 'CHROM', 'ANNOTSV_SCORE', and 'CAUSAL' columns) are processed.
    """
    all_data = []

    # Loop over all files in the directory
    for filename in os.listdir(input_dir):
        file_path = os.path.join(input_dir, filename)

        # Only process files with '.tsv' extension
        if filename.endswith('.tsv'):
            try:
                # Read the file into a DataFrame
                df = pd.read_csv(file_path, sep='\t')

                # Check if the required columns exist
                required_columns = ['CHROM', score, 'CAUSAL']
                if all(col in df.columns for col in required_columns):
                    # Filter out INS for AnnotSV
                    if score == 'ANNOTSV_SCORE' or score == 'ANNOTSV_ACMG_CLASS':
                        if 'TYPE' in df.columns:
                            df = df[df['TYPE'] != 'INS']
                        all_data.append(df)
                    else:
                        all_data.append(df)
                else:
                    print(f"Skipping file {filename}: Missing required columns")
            except Exception as e:
                print(f"Error reading {filename}: {e}")

    # Combine all DataFrames into one large DataFrame
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        return combined_df
    else:
        raise ValueError("No valid files found in the directory.")


def create_dot_plot(input_dir, output_plot, name, score):
    """
    Dotplot with all dots aligned on a single vertical line (scatter).

    Color variants based on the CAUSAL column:
    - 'N': Grey
    - 'Y' or 'Y*': Blue
    - 'y': Orange
    """
    df = read_and_combine_files(input_dir, score)

    # Ensure necessary columns exist
    required_columns = ['CHROM', score, 'CAUSAL']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Input file must contain the following columns: {', '.join(required_columns)}")

    # Create a plot
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    # Masks for different CAUSAL values
    causal_n_mask = df['CAUSAL'] == 'N'
    causal_y_upper_mask = df['CAUSAL'] == 'Y'
    causal_y_middle_mask = df['CAUSAL'] =='Y*'
    causal_y_lower_mask = df['CAUSAL'] == 'y'

    # Add jitter to x-axis for better visualization
    x_base = 1
    df['X_JITTERED'] = x_base + np.random.normal(0, 0.05, size=len(df))

    # Plot non-causal variants (N)
    plt.scatter(
        df['X_JITTERED'][causal_n_mask],
        df[score][causal_n_mask],
        s=10,
        alpha=0.4,
        c='slategray',
        edgecolors='dimgrey',
        label='Non-causal (N)'
    )

    # Plot causal variants Y
    plt.scatter(
        df['X_JITTERED'][causal_y_upper_mask],
        df[score][causal_y_upper_mask],
        s=50,
        alpha=1,
        c='dodgerblue',
        edgecolors='royalblue',
        label='Causal (Y)'
    )

    # Plot causal variants Y*
    plt.scatter(
        df['X_JITTERED'][causal_y_middle_mask],
        df[score][causal_y_middle_mask],
        s=50,
        alpha=1,
        c='hotpink',
        edgecolors='deeppink',
        label='Causal (Y*)'
    )

    # Plot causal variants (y)
    plt.scatter(
        df['X_JITTERED'][causal_y_lower_mask],
        df[score][causal_y_lower_mask],
        s=50,
        alpha=0.8,
        c='orange',
        edgecolors='darkorange',
        label='Causal (y)'
    )

    plt.xlabel("")
    plt.ylabel(score)
    plt.title(f"{score} scores for WGS data")

    if score == 'ANNOTSV_SCORE' :
        plt.ylim(-2, 2)
        custom_ticks = [-1.5, -0.99, -0.9, 0, 0.9, 0.99, 1.5]
        plt.yticks(custom_ticks)
        plt.axhline(0.99, color='r', lw=1, linestyle='--')
        plt.axhline(-0.99, color='r', lw=1, linestyle='--')
        plt.axhline(0.90, color='C1', lw=1, linestyle=':')
        plt.axhline(-0.90, color='C1', lw=1, linestyle=':')
    elif score == 'MAX_PATH_SCORE' or score == 'MAX_OVERLAP_SCORE' :
        plt.ylim(-1, 50)
        plt.axhline(10, color='r', lw=1, linestyle='--')

    # # ALL ACMG CLASSES OF SCORED CAUSAL VARIANTS
    # causalClassSubset = df.query("CAUSAL==['Y','y','Y*']")['ANNOTSV_ACMG_CLASS']
    # causalFrame = causalClassSubset.to_frame(name="ACMG_CLASS")
    # # print(causalFrame)
    #
    # try:
    #     print('Class 5: ' + str(causalFrame['ACMG_CLASS'].value_counts()[5]))
    #     print('Class 4: ' + str(causalFrame['ACMG_CLASS'].value_counts()[4]))
    #     print('Class 3: ' + str(causalFrame['ACMG_CLASS'].value_counts()[3]))
    #     # print('Class 2: ' + str(causalFrame['ACMG_CLASS'].value_counts()[2]))
    #     print('Class 1: ' + str(causalFrame['ACMG_CLASS'].value_counts()[1]))
    # except KeyError:
    #     print('Class not present')

    ## AMOUNT OF CAUSAL VARIANTS: TOTAL AND SCORED
    if score == 'MAX_PATH_SCORE' or score == 'MAX_OVERLAP_SCORE':
        causal_amount = df[df["CAUSAL"].isin(['Y', 'y', 'Y*'])]
        causal_amount_scored = df[df["CAUSAL"].isin(['Y', 'y', 'Y*']) & df[score].notna()]
        causal_above_10 = df[
            df["CAUSAL"].isin(['Y', 'y', 'Y*']) & df[score].notna() & (df[score] > 15)]
        print(f"Scored causal variants with a score above 10: {causal_above_10.shape[0]}")
    else:
        causal_amount = df[df["CAUSAL"].isin(['Y', 'y', 'Y*'])]
        causal_amount_scored = df[df["CAUSAL"].isin(['Y', 'y', 'Y*']) & df["ANNOTSV_SCORE"].notna()]
    print(f"Causal amount total: {causal_amount.shape[0]}")
    print(f"Causal amount Scored: {causal_amount_scored.shape[0]}")

    ## AMOUNT OF NON-CAUSAL VARIANTS: TOTAL AND SCORED
    if score == 'MAX_PATH_SCORE':
        noncausal_amount = df[~df["CAUSAL"].isin(['Y', 'y', 'Y*'])]
        noncausal_amount_scored = df[~df["CAUSAL"].isin(['Y', 'y', 'Y*']) & df[score].notna()]
        noncausal_above_10 = noncausal_amount[noncausal_amount[score] > 14.99]
        print(f"Scored NON causal variants with a score above 14.99: {noncausal_above_10.shape[0]}")
        print(f"NON Causal amount total: {noncausal_amount.shape[0]}")
        print(f"NON Causal amount scored: {noncausal_amount_scored.shape[0]}")
    elif score == 'MAX_OVERLAP_SCORE':
        noncausal_amount = df[~df["CAUSAL"].isin(['Y', 'y', 'Y*']) & (df["TYPE"] != "INS")]
        noncausal_amount_scored = noncausal_amount[noncausal_amount[score].notna()]
        noncausal_above_10 = noncausal_amount[noncausal_amount[score] > 9.99]
        print(f"Scored NON causal variants with a score above 9.99: {noncausal_above_10.shape[0]}")
        print(f"NON Causal amount total: {noncausal_amount.shape[0]}")
        print(f"NON Causal amount scored: {noncausal_amount_scored.shape[0]}")
    else:
        noncausal_amount = df[~df["CAUSAL"].isin(['Y', 'y', 'Y*'])]
        noncausal_amount_scored = df[~df["CAUSAL"].isin(['Y', 'y', 'Y*']) & df[score].notna()]
        noncausal_above_10 = noncausal_amount[noncausal_amount[score] > 0.98]
        print(f"Scored NON causal variants with a score above 0.98: {noncausal_above_10.shape[0]}")
        print(f"NON Causal amount total: {noncausal_amount.shape[0]}")
        print(f"NON Causal amount scored: {noncausal_amount_scored.shape[0]}")


    # Set y-axis limit
    if score == 'ANNOTSV_SCORE':
        plt.ylim(-2, 2)
    elif score == 'ANNOTSV_ACMG_CLASS':
        plt.ylim(0.5, 5.5)
    else:
        plt.ylim(-0.5, 50)

    # Remove x-axis ticks and labels
    plt.xticks([])
    plt.grid(axis='y', linestyle='--', alpha=0.4)

    # Save the plot
    plt.tight_layout()
    plt.legend()
    # plt.savefig(output_plot)
    # plt.show()

    print(f"Dot plot saved as: {output_plot}")



if __name__ == "__main__":

    caller = "LRSR"
    # score = "MAX_OVERLAP_SCORE"
    score = "MAX_PATH_SCORE"
    # score = "ANNOTSV_SCORE"
    # score = "ANNOTSV_ACMG_CLASS"


    inputdir="C:/Users/Tessa.vanderVeer@radboudumc.nl/Documents/Result_analysis/"+caller+"_data/processed_CT"
    output_plot="C:/Users/Tessa.vanderVeer@radboudumc.nl/Documents/Result_analysis/"+caller+"_data/" + caller + "_" + score + "_noChr_dot_plot.png"

    create_dot_plot(inputdir, output_plot, caller, score)
