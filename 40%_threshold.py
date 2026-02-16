import pandas as pd

df = pd.read_csv("protac_neurodegenerative_final.csv")

# Clean degradation column again
df["Percent degradation (%)"] = df["Percent degradation (%)"].astype(str)
df["Percent degradation (%)"] = df["Percent degradation (%)"].str.extract(r'(\d+\.?\d*)')
df["Percent degradation (%)"] = pd.to_numeric(df["Percent degradation (%)"], errors="coerce")

# Use 40% threshold
df["label"] = df["Percent degradation (%)"].apply(
    lambda x: 1 if x >= 40 else 0
)

print(df["label"].value_counts())

df.to_csv("protac_neurodegenerative_labeled_40percent.csv", index=False)
