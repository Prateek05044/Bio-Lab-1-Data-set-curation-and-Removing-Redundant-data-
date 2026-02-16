import pandas as pd

# Load dataset
df = pd.read_csv("protac_neurodegenerative_labeled.csv")

print("Initial dataset size:", df.shape)

# --------------------------------------------------
# 1️⃣ Remove completely identical rows
# --------------------------------------------------
df = df.drop_duplicates()
print("After removing exact duplicate rows:", df.shape)

# --------------------------------------------------
# 2️⃣ Remove duplicate canonical SMILES (same molecule)
# --------------------------------------------------
if "canonical_smiles" in df.columns:
    before = df.shape[0]
    df = df.drop_duplicates(subset=["canonical_smiles"])
    after = df.shape[0]
    print(f"Removed {before - after} duplicate molecules (canonical SMILES)")

# --------------------------------------------------
# 3️⃣ Remove duplicate molecule-target pairs
# --------------------------------------------------
if "canonical_smiles" in df.columns and "Target" in df.columns:
    before = df.shape[0]
    df = df.drop_duplicates(subset=["canonical_smiles", "Target"])
    after = df.shape[0]
    print(f"Removed {before - after} duplicate molecule-target pairs")

print("Final dataset size after redundancy removal:", df.shape)

# --------------------------------------------------
# 4️⃣ Save cleaned dataset
# --------------------------------------------------
df.to_csv("protac_neurodegenerative_cleaned.csv", index=False)

print("✅ Cleaned dataset saved as: protac_neurodegenerative_cleaned.csv")
