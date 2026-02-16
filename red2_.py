import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

# 1. Load the dataset
file_name = 'protac_neurodegenerative_cleaned.csv'
df = pd.read_csv(file_name)

print(f"--- Dataset Integrity Check for: {file_name} ---")

# 2. Redundancy Check (Exact Matches)
# Checking for duplicate structures using Canonical SMILES and InChI Keys
dup_smiles = df.duplicated(subset=['canonical_smiles']).sum()
dup_inchi = df.duplicated(subset=['InChI Key']).sum()

print(f"Duplicate Canonical SMILES found: {dup_smiles}")
print(f"Duplicate InChI Keys found: {dup_inchi}")

# 3. Tanimoto Similarity Test
print("\nCalculating Tanimoto Similarity (this may take a minute)...")

# Convert SMILES to RDKit Molecule objects
mols = [Chem.MolFromSmiles(s) for s in df['canonical_smiles']]

# Generate Morgan Fingerprints (Radius 3 is equivalent to ECFP6)
# We use useChirality=True to distinguish between stereoisomers
fps = [AllChem.GetMorganFingerprintAsBitVect(m, 3, nBits=2048, useChirality=True) for m in mols if m is not None]

# Calculate pairwise Tanimoto Similarity
n = len(fps)
similarities = []
for i in range(n):
    for j in range(i + 1, n):
        sim = DataStructs.FingerprintSimilarity(fps[i], fps[j])
        similarities.append(sim)

# 4. Statistical Summary
avg_sim = np.mean(similarities)
max_sim = np.max(similarities)
median_sim = np.median(similarities)

print(f"Average Tanimoto Similarity: {avg_sim:.4f}")
print(f"Median Tanimoto Similarity: {median_sim:.4f}")
print(f"Maximum Tanimoto Similarity: {max_sim:.4f}")

# 5. Visualization
plt.figure(figsize=(10, 6))
sns.histplot(similarities, bins=50, kde=True, color='teal')
plt.axvline(avg_sim, color='red', linestyle='--', label=f'Mean: {avg_sim:.2f}')
plt.title('Structural Diversity: Distribution of Pairwise Tanimoto Similarities')
plt.xlabel('Tanimoto Similarity (0 = Dissimilar, 1 = Identical)')
plt.ylabel('Number of Molecule Pairs')
plt.legend()

# Save the plot for your report
plt.savefig('tanimoto_report_plot.png', dpi=300)
print("\nSuccess! Plot saved as 'tanimoto_report_plot.png'.")

# 6. Identifying Near-Duplicates (Sim > 0.90)
high_sim_count = sum(1 for s in similarities if s > 0.90)
print(f"Number of highly similar pairs (Sim > 0.90): {high_sim_count}")