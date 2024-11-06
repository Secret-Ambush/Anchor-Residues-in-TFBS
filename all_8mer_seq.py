import itertools

# Function to generate complement of a DNA sequence
def complement(seq):
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement_map[base] for base in seq)

# Generate all possible 8-length combinations of DNA
bases = ['A', 'T', 'C', 'G']
combinations = [''.join(p) for p in itertools.product(bases, repeat=8)]

# Write output to a text file
with open('dna_combinations.txt', 'w') as file:
    file.write("Sequence\tComplement\n")
    for seq in combinations:
        file.write(f"{seq}\t{complement(seq)}\n")

print("DNA combinations and complements have been written to 'dna_combinations.txt'.")