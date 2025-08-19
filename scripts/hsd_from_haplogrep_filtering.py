import pandas as pd
import re
from pathlib import Path

# Load the CSV file
file_path = "~/Downloads/haplogroups_consensus/haplogroups.extended.csv"  # adjust path if needed
df = pd.read_csv(file_path, sep="\t")

# Keep only relevant columns
df = df[["SampleID", "Input_Sample"]]

# Function to extract positions for range calculation
def extract_positions(mutations_str):
    positions = []
    if pd.isna(mutations_str):
        return []
    for mut in mutations_str.split():
        # accept if:
        # 1. ends with canonical base (C,T,A,G)
        # 2. OR contains insertion/deletion (".", "d")
        if mut[-1] in ["C", "T", "A", "G"] or "d" in mut or "." in mut:
            match = re.match(r"(\d+)", mut)
            if match:
                positions.append(int(match.group(1)))
    return positions

# Normalize sample name to group HVS regions together
def normalize_sample(sample_name):
    return re.sub(r"_HVS\d+_consensus", "", sample_name)

df["SampleBase"] = df["SampleID"].apply(normalize_sample)

# Group by SampleBase
grouped = []
for sample, subdf in df.groupby("SampleBase"):
    all_mutations = []
    all_ranges = []
    regions = []
    for _, row in subdf.iterrows():
        # Extract region type
        match_region = re.search(r"(HVS\d+)", row["SampleID"])
        if match_region:
            regions.append(match_region.group(1))

        if pd.notna(row["Input_Sample"]):
            # keep valid mutations (substitutions + indels)
            valid_muts = [m for m in row["Input_Sample"].split()
                          if (m[-1] in ["C", "T", "A", "G"] or "d" in m or "." in m)]
            all_mutations.extend(valid_muts)

            # extract positions for ranges
            positions = extract_positions(" ".join(valid_muts))
            if positions:
                all_ranges.append(f"{min(positions)}-{max(positions)}")

    # Sort regions
    regions_sorted = sorted(regions, key=lambda x: int(x.replace("HVS", "")))
    region_pattern = "_".join(regions_sorted)

    # Construct new sample name
    new_sample_id = f"{sample}_{region_pattern}"

    # Concatenate subranges with ";"
    range_str = ";".join(all_ranges)

    # Construct mutations string
    mutations_str = " ".join(all_mutations)

    grouped.append({
        "SampleID": new_sample_id,
        "Range": range_str,
        "Polymorphisms": mutations_str,
        "Haplogroup": "?"
    })

# Create output dataframe
df_out = pd.DataFrame(grouped, columns=["SampleID", "Range", "Polymorphisms", "Haplogroup"])

# Save to .hsd file
output_path = Path(file_path).with_suffix(".hsd")
df_out.to_csv(output_path, sep="\t", index=False)

print(f"File saved to {output_path}")

