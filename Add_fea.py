import pandas as pd

df = pd.read_csv("protac_neurodegenerative_final.csv")

# Convert degradation column to numeric (just in case)
df["Percent degradation (%)"] = pd.to_numeric(
    df["Percent degradation (%)"], errors="coerce"
)

# Create proper labels
df["label"] = df["Percent degradation (%)"].apply(
    lambda x: 1 if x >= 50 else 0
)

# Check distribution
print("Label distribution:")
print(df["label"].value_counts())

# Save updated dataset
df.to_csv("protac_neurodegenerative_labeled.csv", index=False)

print("Updated labeled dataset saved.")
