import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

# 1. Load the original dataset
input_file = 'protac_neurodegenerative_labeled_40percent.csv'
df = pd.read_csv(input_file)

print(f"Starting with {len(df)} rows.")

# 2. Convert SMILES to Fingerprints (ECFP6 with Chirality)
# We use Radius 3 (ECFP6) and Chirality to be as precise as possible
mols = [Chem.MolFromSmiles(s) for s in df['canonical_smiles']]
fps = [AllChem.GetMorganFingerprintAsBitVect(m, 3, nBits=2048, useChirality=True) for m in mols]

# 3. Identify and filter out redundancies (Similarity = 0.95)
indices_to_keep = []
for i in range(len(fps)):
    is_redundant = False
    # Check current molecule against all previously accepted molecules
    for j in indices_to_keep:
        sim = DataStructs.FingerprintSimilarity(fps[i], fps[j])
        if sim >= 0.95:
            is_redundant = True
            break
    if not is_redundant:
        indices_to_keep.append(i)

# 4. Create the final curated dataset
df_final = df.iloc[indices_to_keep].copy()

# 5. Save to a new Excel file
output_file = 'PROTAC_Tanimoto_Curated.xlsx'
df_final.to_excel(output_file, index=False)

print(f"Removed {len(df) - len(df_final)} redundant structures.")
print(f"Final dataset has {len(df_final)} unique molecules.")
print(f"File saved successfully as: {output_file}")