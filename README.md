# RNA Phage Host Prediction Analysis

This script processes bacterial family data from multiple samples to identify key families based on their relative abundance and consistency across samples. It uses thresholds to progressively filter the data for relevance.

## How It Works

1. **Inner Threshold (`InnerPercentAbove`)**: Filters families within a sample based on relative abundance.
2. **Samples Threshold (`samplesPercentAbove`)**: Filters families across multiple samples based on their prevalence.
3. **Final Threshold (`inAllLinesPercentAbove`)**: Retains families that are consistently found across all results.

## Requirements

- Python 3.x
- `pandas` library
- Excel file with a list of samples
- `.tsv` files containing bacterial family data for each sample

## Usage

1. Update the paths in the script to match your file locations.
2. Set the thresholds (`InnerPercentAbove`, `samplesPercentAbove`, `inAllLinesPercentAbove`) as needed.
3. Run the script to generate a filtered list of families and save the results to a CSV file.

## Example

- A family that appears in at least **13%** of a single sample, **95%** of all samples, and **95%** of the overall results will be included in the final list.

## Output

The script generates a CSV file with the following columns:
- Sample-specific data
- Percentages of families within samples
- Filtered families based on thresholds
- Final consistent families

For more details on the thresholds, see the comments in the script.
