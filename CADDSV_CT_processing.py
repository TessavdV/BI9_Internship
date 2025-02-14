import pandas as pd

"""
Isoleert CT-O en CT-P uit alle CT-scores. 
Isolates CT-O and CT-P score from all CT-scores for all output results. 

Inputfile bevatten Variant, AnnotSV-score en alle gevonden CT-scores (en bijbehorende varianten). 
"""

def calc_max_path(df):
    """
    Calculate the maximum pathogenic score and its associated variant.
    """
    df['MAX_PATH_SCORE'] = None
    df['MAX_PATH_VAR'] = None

    for index, row in df.iterrows():
        vars_value = row['CADDSV_VARS']
        scores_value = row['CADDSV_SCORE']

        if vars_value != "Not Present" and scores_value != "Not Present":
            # print(vars_value)
            caddsv_vars = vars_value.split(",")
            caddsv_scores = []

            # Parse scores, skipping invalid entries (formatting error)
            for score in scores_value.split(","):
                try:
                    caddsv_scores.append(float(score))
                except ValueError:
                    caddsv_scores.append(float('-inf'))

            # Find the most pathogenic score
            if caddsv_scores:
                max_path_score, max_path_var = max(
                    zip(caddsv_scores, caddsv_vars),
                    key=lambda x: x[0]
                )
                if max_path_score != float('-inf'):
                    df.at[index, 'MAX_PATH_SCORE'] = max_path_score
                    df.at[index, 'MAX_PATH_VAR'] = max_path_var

    return df

def calc_max_overlap(df):
    """
    Calculate the variant with the maximum overlap and its associated score, min. 10% reciprocal overlap.
    """
    df['MAX_OVERLAP_SCORE'] = None
    df['MAX_OVERLAP_VAR'] = None

    for index, row in df.iterrows():
        start_pos = row['START']
        end_pos = row['END']
        vars_value = row['CADDSV_VARS']
        scores_value = row['CADDSV_SCORE']

        if vars_value != "Not Present" and scores_value != "Not Present":
            caddsv_vars = vars_value.split(",")
            caddsv_scores = []

            # Parse scores, skipping invalid entries
            for score in scores_value.split(","):
                try:
                    caddsv_scores.append(float(score))
                except ValueError:
                    caddsv_scores.append(float('-inf'))

            # Calculate overlaps
            valid_overlaps = []
            valid_scores = []
            valid_vars = []

            for i, var in enumerate(caddsv_vars):
                try:
                    var_start, var_end = map(int, var.split(":")[1].split("-"))
                    overlap = max(0, min(end_pos, var_end) - max(start_pos, var_start))

                    # Calculate the reciprocal overlaps
                    variant_len = end_pos - start_pos
                    caddsv_len = var_end - var_start

                    reciprocal_overlap_variant = overlap / variant_len if variant_len > 0 else 0
                    reciprocal_overlap_caddsv = overlap / caddsv_len if caddsv_len > 0 else 0

                    if reciprocal_overlap_variant >= 0.1 and reciprocal_overlap_caddsv >= 0.1:
                        valid_overlaps.append(overlap)
                        valid_scores.append(caddsv_scores[i])
                        valid_vars.append(var)

                except Exception:
                    continue

            # Find maximum overlap among valid overlaps
            if valid_overlaps:
                max_overlap_idx = valid_overlaps.index(max(valid_overlaps))
                max_overlap_score = valid_scores[max_overlap_idx]
                df.at[index, 'MAX_OVERLAP_SCORE'] = max_overlap_score
                df.at[index, 'MAX_OVERLAP_VAR'] = valid_vars[max_overlap_idx]

    return df


def process_file(input_file, output_file):
    """
    Main processing function to read the input file, calculate scores, and save the output.
    """
    # Load the file into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t')

    # Ensure necessary columns exist
    required_columns = ['CHROM', 'START', 'END', 'CADDSV_VARS', 'CADDSV_SCORE']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Input file must contain the following columns: {', '.join(required_columns)}")

    # Calculate scores
    df = calc_max_path(df)
    df = calc_max_overlap(df)

    # Remove unnecessary columns
    df.drop(columns=['CADDSV_VARS', 'CADDSV_SCORE'], inplace=True)

    # Save to output file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Processing complete. Updated file saved as '{output_file}'.")


if __name__ == "__main__":
    samplename = "P50-G6"
    input_file = "FILEPATH_HERE" + samplename +  "_fullCADDSV_results.tsv"
    output_file = "FILEPATH_HERE" + samplename + "_CADDSV_CTCPandFT.tsv"

    process_file(input_file, output_file)


