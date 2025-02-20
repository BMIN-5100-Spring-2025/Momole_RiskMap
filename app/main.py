import pandas as pd
import os
import json

# Define input and output folders
INPUT_FOLDER = "input"
OUTPUT_FOLDER = "output"

# Ensure output directory exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

def find_latest_file(directory, extension):
    """Find the latest file in a directory with the given extension."""
    files = [f for f in os.listdir(directory) if f.endswith(extension)]
    if not files:
        return None
    files.sort(key=lambda x: os.path.getmtime(os.path.join(directory, x)), reverse=True)
    return os.path.join(directory, files[0])

def calculate_prs(gwas_path, genotype_path):
    """Computes PRS from GWAS summary and genotype data."""
    print(f"Loading GWAS file: {gwas_path}")
    gwas_df = pd.read_csv(gwas_path, sep="\t")

    # Check for required columns
    required_cols = {"SNP", "Effect_Allele", "Beta"}
    if not required_cols.issubset(gwas_df.columns):
        raise ValueError(f"Missing columns in GWAS file: {required_cols - set(gwas_df.columns)}")

    print(f"Loading Genotype file: {genotype_path}")
    genotype_df = pd.read_csv(genotype_path, sep="\t")

    # Ensure SNPs match in both files
    common_snps = set(gwas_df["SNP"]).intersection(set(genotype_df["SNP"]))
    if len(common_snps) == 0:
        raise ValueError("No matching SNPs between GWAS and genotype data.")

    # Merge DataFrames on SNP
    merged_df = genotype_df.merge(gwas_df, on="SNP")

    # PRS Calculation: Weighted sum of effect alleles
    merged_df["PRS_Component"] = merged_df["Genotype"] * merged_df["Beta"]
    prs_score = merged_df.groupby("Individual_ID")["PRS_Component"].sum().reset_index()
    prs_score.columns = ["Individual_ID", "PRS_Score"]

    return prs_score

def main():
    """Main function to process PRS calculation."""
    # Find latest GWAS and Genotype files
    gwas_file = find_latest_file(INPUT_FOLDER, ".tsv")
    genotype_file = find_latest_file(INPUT_FOLDER, ".tsv")

    if not gwas_file or not genotype_file:
        print("❌ Missing GWAS or Genotype file in 'input/' folder.")
        return

    try:
        prs_results = calculate_prs(gwas_file, genotype_file)

        # Save PRS results
        prs_csv_path = os.path.join(OUTPUT_FOLDER, "prs_results.csv")
        prs_json_path = os.path.join(OUTPUT_FOLDER, "prs_results.json")

        prs_results.to_csv(prs_csv_path, index=False)
        prs_results.to_json(prs_json_path, orient="records")

        print(f"✅ PRS calculation complete! Results saved in:\n📁 {prs_csv_path}\n📁 {prs_json_path}")

    except Exception as e:
        print(f"❌ Error: {str(e)}")

if __name__ == "__main__":
    main()