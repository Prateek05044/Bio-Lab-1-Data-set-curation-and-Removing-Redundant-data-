import pandas as pd
from rdkit import Chem

# --------------------------------------------------
# 1. LOAD FULL PROTAC-DB
# --------------------------------------------------
df = pd.read_csv("protac_db.csv")
print("Initial PROTAC-DB size:", df.shape)

# --------------------------------------------------
# 2. DEFINE NEURODEGENERATIVE TARGET SET
# --------------------------------------------------
neuro_targets = [
    # Core neurodegenerative disease proteins
    "MAPT", "Tau",
    "SNCA", "alpha-synuclein",
    "HTT", "huntingtin",
    "TDP-43", "TARDBP",
    "APP", "amyloid",
    "SOD1",

    # Protein quality control & neurodegeneration pathways
    "PARK", "parkin",
    "LRRK2",
    "UBQLN",
    "VCP",
    "PSEN",
    "HSP", "chaperone",
    "proteasome",
    "HDAC"
]

pattern = "|".join(neuro_targets)

# --------------------------------------------------
# 3. FILTER FOR NEURODEGENERATIVE TARGETS
# --------------------------------------------------
df_neuro = df[df["Target"].str.contains(
    pattern, case=False, na=False
)]

print("After neurodegenerative target filtering:", df_neuro.shape)

# --------------------------------------------------
# 4. DROP MISSING SMILES
# --------------------------------------------------
df_neuro = df_neuro.dropna(subset=["Smiles"])

# --------------------------------------------------
# 5. VALIDATE SMILES (RDKIT)
# --------------------------------------------------
df_neuro["mol"] = df_neuro["Smiles"].apply(Chem.MolFromSmiles)
df_neuro = df_neuro[df_neuro["mol"].notnull()].reset_index(drop=True)

print("After SMILES validation:", df_neuro.shape)

# --------------------------------------------------
# 6. CANONICALIZE SMILES
# --------------------------------------------------
df_neuro["canonical_smiles"] = df_neuro["mol"].apply(
    lambda mol: Chem.MolToSmiles(mol, canonical=True)
)

# --------------------------------------------------
# 7. REMOVE REDUNDANT PROTACS
# --------------------------------------------------
before = df_neuro.shape[0]
df_neuro = df_neuro.drop_duplicates(subset=["canonical_smiles"])
after = df_neuro.shape[0]

print(f"Removed {before - after} duplicate PROTACs")
print("After deduplication:", df_neuro.shape)

# --------------------------------------------------
# 8. OPTIONAL LABEL (SCREENING SETUP)
# --------------------------------------------------
if "degradation_percent" in df_neuro.columns:
    df_neuro["label"] = df_neuro["degradation_percent"].apply(
        lambda x: 1 if x >= 50 else 0
    )
else:
    df_neuro["label"] = 0  # unknown / screening

# --------------------------------------------------
# 9. CLEANUP
# --------------------------------------------------
df_neuro = df_neuro.drop(columns=["mol"])

# --------------------------------------------------
# 10. SAVE FINAL DATASET
# --------------------------------------------------
df_neuro.to_csv(
    "protac_neurodegenerative_final.csv",
    index=False
)

print("✅ FINAL NEURODEGENERATIVE DATASET SAVED")
print("Final dataset size:", df_neuro.shape)
