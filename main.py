import pandas as pd
import numpy as np

#1. LOCATING THE DATA & METADATA
# Raw bioinformatics files have a lot of text notes.
# I need to find the actual table start and the lines with patient labels.
filename = "GSE19188_series_matrix.txt"
print("Checking the file")

with open(filename, 'r') as file_data:
    lines = file_data.readlines()

start_row = 0
metadata_lines = []

for index, text in enumerate(lines):
    # I grab lines starting with '!' because the 'tumor' labels are hidden there.
    if text.startswith('!Sample_'):
        metadata_lines.append(text.lower())
    # This marker tells me exactly where the table of numbers begins.
    if "!series_matrix_table_begin" in text:
        start_row = index
        break

print("The real data starts at line:", start_row)

#2. LOADING THE TABLE 
# I skip the intro lines and tell Python the data is separated by tabs.
raw_table = pd.read_csv(filename, skiprows=start_row + 1, sep='\t', comment='!')

# Cleaning: Remove empty rows and flip it (transpose) 
# so each row is a Patient and each column is a Gene.
raw_table = raw_table.dropna(how='all')
raw_table = raw_table.set_index('ID_REF')
clean_table = raw_table.T 

#3. CONVERTING TO NUMBERS
# Sometimes data is read as text. This forces it to be decimal numbers.
clean_table = clean_table.apply(pd.to_numeric, errors='coerce')

# If any values are missing, I fill them with 0.
# I am using 0 because calculating 54,000 medians was taking too long 
# and making the computer hang.
clean_table = clean_table.fillna(0)

#4. LOG2 TRANSFORMATION 
# Gene numbers are too big and uneven. Log2 squashes them to a 0-16 scale.
# I add 1 to avoid math errors with zero values.
final_data = np.log2(clean_table + 1)

#5. CREATING GROUPS (TUMOR VS NORMAL)
# I create a list of 'Normal' labels first for all patients.
num_patients = len(raw_table.columns)
group_labels = ["Normal"] * num_patients 

# I search every note I saved in Step 1. If 'tumor' appears in a patient's 
# description, I change their label to 'Tumor'.
for line in metadata_lines:
    # I split the line by tabs and skip the first item because it's just the row title.
    parts = line.replace('"', '').split('\t')[1:]
    for i in range(len(parts)):
        if i < len(group_labels):
            if 'tumor' in parts[i]:
                group_labels[i] = "Tumor"

print("Found labels:", set(group_labels))
final_data['Target'] = group_labels

#6. FIXING THE POWER BI LIMIT 
# Power BI crashes with 54,000 columns. I 'melt' them into one long list.
print("Melting data for Power BI...")
output_for_bi = final_data.melt(id_vars=['Target'], var_name='Gene_ID', value_name='Expression')

#7. SAVING THE PROJECT FOR FURTHER ANALYSIS
output_for_bi.to_csv("Lung_Cancer_Analysis_Final.csv", index=False)
print("Project ready for analysis!")
