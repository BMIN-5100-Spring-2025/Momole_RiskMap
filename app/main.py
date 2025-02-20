import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import geopandas as gpd

# Set input and output folders
INPUT_FOLDER = "input"
OUTPUT_FOLDER = "output"

# Ensure output directory exists
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# Expected column mappings (now using Gene Set instead of SNP)
COLUMN_MAPPINGS = {
    "Gene_Set": ["Gene Set", "Geneset", "Pathway", "Curated Set"],
    "Beta": ["Beta", "Effect_Size"],
}

def find_latest_file(directory, extensions):
    """Find the latest file in a directory with the given extension."""
    files = [f for f in os.listdir(directory) if f.endswith(tuple(extensions))]
    if not files:
        return None
    files.sort(key=lambda x: os.path.getmtime(os.path.join(directory, x)), reverse=True)
    return os.path.join(directory, files[0])

def standardize_columns(df):
    """Rename GWAS columns to match expected names."""
    rename_map = {}
    for standard_name, possible_names in COLUMN_MAPPINGS.items():
        for col in possible_names:
            if col in df.columns:
                rename_map[col] = standard_name
                break
    df.rename(columns=rename_map, inplace=True)
    return df

def detect_separator(file_path):
    """Detect if file is CSV (comma-separated) or TSV (tab-separated)."""
    with open(file_path, "r", encoding="utf-8") as f:
        first_line = f.readline()
        return "," if "," in first_line else "\t"

def calculate_gene_based_prs(gwas_path):
    """Computes PRS from gene-based GWAS summary data."""
    print(f"Loading GWAS file: {gwas_path}")

    # Detect file format
    separator = detect_separator(gwas_path)
    gwas_df = pd.read_csv(gwas_path, sep=separator)

    # Standardize column names
    gwas_df = standardize_columns(gwas_df)

    # Check if required columns exist after renaming
    required_cols = {"Gene_Set", "Beta"}
    if not required_cols.issubset(gwas_df.columns):
        raise ValueError(f"Missing required columns in GWAS file after renaming: {required_cols - set(gwas_df.columns)}")

    # Compute PRS Score for each gene set (sum of beta values)
    prs_score = gwas_df.groupby("Gene_Set")["Beta"].sum().reset_index()
    prs_score.columns = ["Gene_Set", "PRS_Score"]

    return prs_score

def plot_prs_distribution(prs_results):
    """Generates a bar chart for PRS scores."""
    plt.figure(figsize=(10, 5))
    prs_results.plot(kind="bar", x="Gene_Set", y="PRS_Score", legend=False, color="skyblue")
    plt.xlabel("Gene Set")
    plt.ylabel("PRS Score")
    plt.title("Polygenic Risk Score (PRS) Distribution by Gene Set")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    # Save the plot
    plot_path = os.path.join(OUTPUT_FOLDER, "prs_distribution.png")
    plt.savefig(plot_path)
    plt.close()
    print(f"📊 PRS distribution plot saved as {plot_path}")

def generate_risk_map(prs_results):
    """Generates a simulated geographic risk map for PRS scores."""
    world = gpd.read_file("https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip")

    # Simulated PRS mapping to random countries
    country_mapping = {
        "Curated_g": "USA",
        "Other_gene_set": "Canada"
    }

    # Create dataset with country-level PRS scores
    risk_data = pd.DataFrame(list(country_mapping.items()), columns=["Gene_Set", "Country"])
    risk_data = risk_data.merge(prs_results, on="Gene_Set", how="left")

    # Merge with world map
    world = world.merge(risk_data, left_on="ADMIN", right_on="Country", how="left")


    # Plot the risk map
    fig, ax = plt.subplots(figsize=(12, 6))
    world.plot(column="PRS_Score", cmap="Reds", linewidth=0.8, edgecolor="black", legend=True, ax=ax)
    plt.title("Global Distribution of Polygenic Risk Score")
    plt.savefig(os.path.join(OUTPUT_FOLDER, "prs_risk_map.png"))
    plt.close()
    print(f"🗺️ PRS risk map saved as {OUTPUT_FOLDER}/prs_risk_map.png")

def main():
    """Main function to process PRS calculation for gene-based GWAS data and generate visualizations."""
    gwas_file = find_latest_file(INPUT_FOLDER, [".tsv", ".csv"])

    if not gwas_file:
        print("❌ Missing GWAS file in 'input/' folder.")
        return

    try:
        prs_results = calculate_gene_based_prs(gwas_file)

        # Save PRS results
        prs_csv_path = os.path.join(OUTPUT_FOLDER, "gene_prs_results.csv")
        prs_json_path = os.path.join(OUTPUT_FOLDER, "gene_prs_results.json")

        prs_results.to_csv(prs_csv_path, index=False)
        prs_results.to_json(prs_json_path, orient="records")

        print(f"✅ Gene-based PRS calculation complete! Results saved in:\n📁 {prs_csv_path}\n📁 {prs_json_path}")

        # Generate visual outputs
        plot_prs_distribution(prs_results)
        generate_risk_map(prs_results)

    except Exception as e:
        print(f"❌ Error: {str(e)}")

if __name__ == "__main__":
    main()