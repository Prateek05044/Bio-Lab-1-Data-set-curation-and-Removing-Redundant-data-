import pandas as pd
import numpy as np

df = pd.read_csv("protac_neurodegenerative_final.csv")

# Convert degradation column to string first
df["Percent degradation (%)"] = df["Percent degradation (%)"].astype(str)

# Extract numeric part (remove fractions and text)
df["Percent degradation (%)"] = df["Percent degradation (%)"].str.extract(r'(\d+\.?\d*)')

# Convert to numeric
df["Percent degradation (%)"] = pd.to_numeric(
    df["Percent degradation (%)"], errors="coerce"
)

print(df["Percent degradation (%)"].describe())

print("\nValue counts (binned):")
print(
    pd.cut(
        df["Percent degradation (%)"],
        bins=[0,10,20,30,40,50,60,70,80,90,100]
    ).value_counts()
)
