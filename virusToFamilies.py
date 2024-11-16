# Import necessary libraries
import pandas as pd
from collections import Counter

# Set thresholds for filtering families based on their relative abundance and consistency
InnerPercentAbove = 13  # Percentage threshold for filtering families within a sample
samplesPercentAbove = 95  # Percentage threshold for filtering families across multiple samples
inAllLinesPercentAbove = 95  # Final threshold for filtering families consistently across all rows

# Load the sample list from an Excel file
samples = pd.read_excel("C:/Users/UriGNB02/Documents/host_phages/inputs/ANI90_0015689.xlsx")['Sample'].tolist()

# Initialize dictionaries to store data about families and their counts
num_families = {}
family_counts = {}

# Process each sample to extract bacterial family data
for sample in samples:
    tsv_file = f'C:/Users/UriGNB02/Documents/host_phages/phylodists/{sample}.contigLin.assembled.tsv'
    family_names = []
    with open(tsv_file, "r") as tsvfile:
        for line in tsvfile:
            # Filter for lines that include 'Bacteria' and contain sufficient data
            if 'Bacteria;' in line and line.count(';') > 3:
                if line.count(';') == 4:
                    family_field = line.strip().split(';')[4].split(',')[0]
                else:
                    family_field = line.strip().split(';')[4].split(';')[0]
                family_name = family_field.split('  ')[-1].strip().replace('"', '')

                # Include only classified families with a score greater than 0.5
                if family_name != "unclassified" and float(line.split(',')[-1]) > 0.5:
                    family_names.append(family_name)

    # Record the number of families and their counts for each sample
    num_families[sample] = len(family_names)
    family_counts[sample] = Counter(family_names)

# Create a DataFrame to store the results
results = pd.DataFrame({
    'Sample': samples,
    'num_families': [num_families[sample] for sample in samples],
    'family_counts': [family_counts[sample] for sample in samples]
})

# Calculate the percentage of each family within each sample
results['family_percentages'] = results.apply(
    lambda row: ', '.join([f"{f}: {count / row['num_families'] * 100:.2f}%" for f, count in row['family_counts'].items()]),
    axis=1
)

# Filter families based on the InnerPercentAbove threshold
results[f'families_above_{InnerPercentAbove}%_inner'] = results.apply(
    lambda row: ', '.join([
        f"{f}: {p}" for f, p in [
            family_perc.split(': ') for family_perc in row['family_percentages'].split(', ')
            if len(family_perc.split(': ')) > 1 and float(family_perc.split(': ')[1][:-1]) > InnerPercentAbove
        ]
    ]),
    axis=1
)

# Count the occurrences of each family across samples
family_samples = {}
for sample, counts in family_counts.items():
    for family in counts.keys():
        if family not in family_samples:
            family_samples[family] = 1
        else:
            family_samples[family] += 1

# Calculate the percentage of samples each family is present in
family_samples_count = results.apply(
    lambda row: [
        f"{f}: {(family_samples[f] / len(samples)) * 100:.2f}%"
        for f, count in row['family_counts'].items()
        if f in [family_perc.split(': ')[0] for family_perc in row[f'families_above_{InnerPercentAbove}%_inner'].split(', ')]
    ],
    axis=1
)
results['family_in_samples_count'] = family_samples_count.apply(lambda x: ', '.join(x))

# Filter families based on the samplesPercentAbove threshold
results[f'filtered_families_in_samples_above_{samplesPercentAbove}%'] = results['family_in_samples_count'].apply(
    lambda x: ', '.join([
        f for f in x.split(', ')
        if len(f.split(': ')) > 1 and float(f.split(': ')[1][:-1]) > samplesPercentAbove
    ])
)

# Perform final filtering based on the inAllLinesPercentAbove threshold
filtered_families_counts = Counter()
for families in results[f'filtered_families_in_samples_above_{samplesPercentAbove}%']:
    for family in families.split(', '):
        family_name = family.split(': ')[0]
        filtered_families_counts[family_name] += 1

total_rows = len(results)
filtered_families_above_threshold = [
    family for family, count in filtered_families_counts.items()
    if count / total_rows >= inAllLinesPercentAbove / 100
]

# Add final filtered families to the results
results[f'filtered_families_All_lines_above{inAllLinesPercentAbove}%'] = ', '.join(filtered_families_above_threshold)

# Save the results to a CSV file
results.to_csv("C:/Users/UriGNB02/Documents/host_phages/outputs/ANI90_0015689_95_13.csv", index=False)
