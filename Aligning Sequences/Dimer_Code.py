import pandas as pd
from weblogo import *

def read_sequences_from_file(file_path, col_index):
    with open(file_path, 'r') as file:
        next(file)
        sequences = []
        for line in file:
            columns = line.split()
            sequences.append(columns[col_index])
    return sequences

def ref_seq(col_index):
    input_file_path = f'/Users/bristi/Desktop/Design Project/Working-with-TF/TFs UniProbe Data/Cebpa/Cebpa_8mers_top_enrichment.txt'
    with open(input_file_path, 'r') as file: 
        next(file)
        for line in file:
            columns = line.split()
            sequence = (columns[col_index].replace('.',''))
            return sequence   

def generate_sorted_substrings(s):
    # Generate all possible substrings
    substrings = [s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)]
    
    # Sort the substrings first by reverse length, then alphabetically
    substrings_sorted = sorted(substrings, key=lambda x: (-len(x), x))
    
    # Remove duplicates without using a set
    unique_substrings = []
    for substring in substrings_sorted:
        if substring not in unique_substrings:
            unique_substrings.append(substring)
    
    return unique_substrings

def align_sequences(sequences, reference_sequence):
    ref_index = 0
    substrings = generate_sorted_substrings(reference_sequence)

    columns = [i for i in range(-7, 15)]  
    df = pd.DataFrame(columns=columns)

    for sequence in sequences:
        aligned = False
        if reference_sequence in sequence:
            idx = sequence.find(reference_sequence)
            aligned = True
            entire = True
        else:
            prev_best_number_of_matches = 0
            best_substring_index = -1
                    
            for i in range(0, len(substrings)):
                total_number_of_matches = 0

                if substrings[i] in sequence:                    
                    currect_substring = substrings[i]
                    idx = sequence.find(currect_substring)
                    substring_index_in_seed = reference_sequence.find(currect_substring)
                    offset = substring_index_in_seed - idx
                    
                    if offset < 0:
                        aligned_seed = reference_sequence[:len(reference_sequence) - abs(offset)]
                        aligned_sequence = sequence[abs(offset):]
                        
                    else:
                        aligned_seed = reference_sequence[abs(offset):]
                        aligned_sequence = sequence[:len(sequence) - abs(offset)]
                    
                    aligned_seed_substrings = generate_sorted_substrings(aligned_seed)
                    for substring in aligned_seed_substrings:
                        if substring in aligned_sequence:
                            total_number_of_matches += 1
                    
                    if prev_best_number_of_matches < total_number_of_matches:
                        best_substring_index = i
                        prev_best_number_of_matches = total_number_of_matches
                        
            
            best_substring = substrings[best_substring_index]
            idx = sequence.find(best_substring)
            substring_index_in_seed = reference_sequence.find(best_substring)
            entire = False
            aligned = True       


        if aligned:
            if entire:
                offset = ref_index - idx
            else:
                offset = substring_index_in_seed - idx
            indices = range(offset, offset + len(sequence))
            df = pd.concat([df, pd.DataFrame([list(sequence)], columns=indices)], ignore_index=True)

        
    return df

def sequence_logo_generator(file_path, col_index):
    fin = open(file_path)
    seqs = read_seq_data(fin)
    logodata = LogoData.from_seqs(seqs)
    logooptions = LogoOptions()
    logooptions.title = "A Logo Title"
    logooptions.color_scheme = std_color_schemes['classic']
    logoformat = LogoFormat(logodata, logooptions)
    png_bytes = png_print_formatter(logodata, logoformat)

    with open(f'DimerBinding/SeqLogos/cebpb_{col_index}.png', 'wb') as fout:
        fout.write(png_bytes)

# running for both complement and regular
for col_index in [0,1]:
    input_file_path = f'/Users/bristi/Desktop/Design Project/Working-with-TF/TFs UniProbe Data/Cebpa/Cebpa_8mers_11111111.txt'
    output_file_path = f'DimerBinding/aligned_Cebpa_col_{col_index}.txt'
    
    sequences = read_sequences_from_file(input_file_path, col_index)

    reference_sequence = ref_seq(col_index)
    print(reference_sequence, ' ', col_index)
    if not reference_sequence:
        raise ValueError(f"No sequence containing '{reference_sequence}' was found in the file.")

    df_aligned = align_sequences(sequences, reference_sequence)

    columns_to_extract = [col for col in range(0, 8) if col in df_aligned.columns]
    extracted_data = df_aligned[columns_to_extract]
    extracted_data = extracted_data.fillna('-')

    formatted_data = extracted_data.apply(lambda row: ''.join(row.astype(str)), axis=1)

    with open(output_file_path, 'w') as file:
        file.write(reference_sequence + '\n')
        for line in formatted_data:
            file.write(line + '\n')
            
    sequence_logo_generator(output_file_path, col_index)