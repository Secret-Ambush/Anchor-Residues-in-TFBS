import random

input_file = '/Users/bristi/Desktop/Design Project/Working-with-TF/SELEX/Dataset/SRR5340731.fasta'
output_file = '/Users/bristi/Desktop/Design Project/Working-with-TF/SELEX/Extracted_FASTA_Data.txt'

with open(input_file, 'r') as file:
    lines = file.readlines()

num_lines_to_extract = int(len(lines) * 0.0001)
extracted_lines = random.sample(lines, num_lines_to_extract)

with open(output_file, 'w') as file:
    file.writelines(extracted_lines)

print(f'Extracted {num_lines_to_extract} lines to {output_file}')