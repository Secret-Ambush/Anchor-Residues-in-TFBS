import pandas as pd
from weblogo import *

def read_sequences_from_file(file_path, col_index):
    with open(file_path, 'r') as file:
        next(file)
        sequences = []
        scores = []
        for line in file:
            columns = line.split()
            if col_index == 0:
                sequences.append(columns[col_index])
            else:
                sequences.append(columns[col_index][::-1])
            scores.append(float(columns[2]))
    return sequences, scores

def ref_seq(col_index):
    input_file_path = f'/Users/bristi/Desktop/Design Project/Working-with-TF/TFs UniProbe Data/Cebpa/Cebpa_8mers_top_enrichment.txt'
    with open(input_file_path, 'r') as file: 
        next(file)
        for line in file:
            columns = line.split()
            if col_index == 0:
                sequence = (columns[col_index].replace('.',''))
            else:
                sequence = (columns[col_index].replace('.',''))
                sequence = sequence[::-1]
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
                    current_substring = substrings[i]
                    idx = sequence.find(current_substring)
                    substring_index_in_seed = reference_sequence.find(current_substring)
                    offset = substring_index_in_seed - idx
                    
                    if offset < 0:
                        aligned_seed = reference_sequence[:len(reference_sequence) - abs(offset)]
                        aligned_sequence = sequence[abs(offset):]
                        
                    else:
                        aligned_seed = reference_sequence[abs(offset):]
                        aligned_sequence = sequence[:len(sequence) - abs(offset)]
                    
                    # aligned_seed_substrings = generate_sorted_substrings(aligned_seed)
                    
                    # for substring in aligned_seed_substrings:
                    #     if substring in aligned_sequence:
                    #         total_number_of_matches += 1
                    
                    str1 = aligned_seed
                    str2 = aligned_sequence
                    min_length = min(len(str1), len(str2))
                    
                    total_number_of_matches = sum(1 for i in range(min_length) if str1[i] == str2[i])
                    
                    
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

def sequence_logo_generator(TF, file_path, col_index):
    fin = open(file_path)
    seqs = read_seq_data(fin)
    logodata = LogoData.from_seqs(seqs)
    logooptions = LogoOptions()
    logooptions.title = "A Logo Title"
    logooptions.color_scheme = std_color_schemes['classic']
    logooptions.yaxis_scale = 1.0
    logoformat = LogoFormat(logodata, logooptions)
    png_bytes = png_print_formatter(logodata, logoformat)

    with open(f'/Users/bristi/Desktop/Design Project/Working-with-TF/DimerBinding/SeqLogos/recent_{TF}_{col_index}.png', 'wb') as fout:
        fout.write(png_bytes)

# running for both complement and regular
for TF in ['Cebpb']:
    for col_index in [0]:
        input_file_path = f'/Users/bristi/Desktop/Design Project/Working-with-TF/TFs UniProbe Data/{TF}/new.txt'
        output_file_path = f'/Users/bristi/Desktop/Design Project/Working-with-TF/DimerBinding/new_{TF}_col_{col_index}.txt'
        
        sequences, scores = read_sequences_from_file(input_file_path,col_index)
        
        max_score = max(scores)
        max_index = scores.index(max_score)
        reference_sequence = sequences[max_index]
        
        if col_index == 1:
            reference_sequence = reference_sequence[::-1]               
            
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
                
        sequence_logo_generator(TF, output_file_path, col_index)