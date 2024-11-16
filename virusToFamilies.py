import pandas as pd
from collections import Counter

InnerPercentAbove = 13
samplesPercentAbove = 95
inAllLinesPercentAbove = 95

samples = pd.read_excel("C:/Users/UriGNB02/Documents/host_phages/inputs/ANI90_0015689.xlsx")['Sample'].tolist()
num_families = {}
family_counts = {}

for sample in samples:
    tsv_file = f'C:/Users/UriGNB02/Documents/host_phages/phylodists/{sample}.contigLin.assembled.tsv'
    family_names = []
    with open(tsv_file, "r") as tsvfile:
        for line in tsvfile:
            if 'Bacteria;' in line and line.count(';') > 3:
                if line.count(';') == 4:
                    family_field = line.strip().split(';')[4].split(',')[0]
                else:
                    family_field = line.strip().split(';')[4].split(';')[0]
                family_name = family_field.split('  ')[-1].strip().replace('"', '')

                if family_name != "unclassified" and float(line.split(',')[-1]) > 0.5:
                    family_names.append(family_name)

    num_families[sample] = len(family_names)
    family_counts[sample] = Counter(family_names)

results = pd.DataFrame({'Sample': samples, 'num_families': [num_families[sample] for sample in samples], 'family_counts': [family_counts[sample] for sample in samples]})
results['family_percentages'] = results.apply(lambda row: ', '.join([f"{f}: {count / row['num_families'] * 100:.2f}%" for f, count in row['family_counts'].items()]), axis=1)
results[f'families_above_{InnerPercentAbove}%_inner'] = results.apply(lambda row: ', '.join([f"{f}: {p}" for f, p in [family_perc.split(': ') for family_perc in row['family_percentages'].split(', ') if len(family_perc.split(': ')) > 1 and float(family_perc.split(': ')[1][:-1]) > InnerPercentAbove]]), axis=1)

family_samples = {}
for sample, counts in family_counts.items():
    for family in counts.keys():
        if family not in family_samples:
            family_samples[family] = 1
        else:
            family_samples[family] += 1

family_samples_count = results.apply(lambda row: [f"{f}: {(family_samples[f] / len(samples)) * 100:.2f}%" for f, count in row['family_counts'].items() if f in [family_perc.split(': ')[0] for family_perc in row[f'families_above_{InnerPercentAbove}%_inner'].split(', ')]], axis=1)
results['family_in_samples_count'] = family_samples_count.apply(lambda x: ', '.join(x))
results[f'filtered_families_in_samples_above_{samplesPercentAbove}%'] = results['family_in_samples_count'].apply(lambda x: ', '.join([f for f in x.split(', ') if len(f.split(': ')) > 1 and float(f.split(': ')[1][:-1]) > samplesPercentAbove]))

filtered_families_counts = Counter()
for families in results[f'filtered_families_in_samples_above_{samplesPercentAbove}%']:
    for family in families.split(', '):
        family_name = family.split(': ')[0]
        filtered_families_counts[family_name] += 1

total_rows = len(results)
filtered_families_above_threshold = [family for family, count in filtered_families_counts.items() if count / total_rows >= inAllLinesPercentAbove / 100]

results[f'filtered_families_All_lines_above{inAllLinesPercentAbove}%'] = ', '.join(filtered_families_above_threshold)
results.to_csv("C:/Users/UriGNB02/Documents/host_phages/outputs/ANI90_0015689_95_13.csv", index=False)



