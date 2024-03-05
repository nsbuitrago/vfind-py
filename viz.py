import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load the CSV file
df = pd.read_csv('30-982598147/17L_variants.csv')
# drop any sequences that contain * in peptide sequence
df = df[~df['peptide'].str.contains('\*')]
# combine the same peptide sequences and sum the count
df = df.groupby('peptide').agg({'count': 'sum'}).reset_index()

breakpoint()

# remove anything with less than 10 counts
df = df[df['count'] >= 10]
breakpoint()

# Define all possible amino acids (adjust as needed)
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

# Initialize a dictionary to hold the frequency data
frequency_matrix = {aa: [0]*len(df['peptide'].iloc[0]) for aa in amino_acids}

# Calculate the frequency of each amino acid at each position
for index, row in df.iterrows():
    sequence = row['peptide']
    read_count = row['count']
    for position, amino_acid in enumerate(sequence):
        if amino_acid in frequency_matrix:
            frequency_matrix[amino_acid][position] += read_count

# Convert the frequency matrix to a DataFrame for easier plotting
frequency_df = pd.DataFrame(frequency_matrix)

# Normalize the frequencies by the total read counts to get relative frequencies
frequency_df = frequency_df.div(frequency_df.sum(axis=1), axis=0)

# Plotting the heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(frequency_df.transpose(), cmap='Blues', annot=False)
plt.title('Amino Acid Frequency by Position')
plt.xlabel('Position')
plt.ylabel('Amino Acid')
plt.show()
