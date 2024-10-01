def read_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary with sequence headers as keys
    and sequences as values.
    
    :param file_path: Path to the FASTA file.
    :return: A dictionary with headers and sequences.
    """
    fasta_dict = {}
    header = None
    sequence = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:  # Save the previous header and sequence
                    fasta_dict[header] = ''.join(sequence)
                header = line[1:]  # Remove the '>' from the header
                sequence = []  # Reset the sequence list for the new header
            else:
                sequence.append(line)
        if header:  # Don't forget to save the last sequence
            fasta_dict[header] = ''.join(sequence)
    
    return fasta_dict

def find_orfs(sequence, min_length=20):
    """
    Finds all ORFs in the DNA sequence that are greater than or equal to the
    specified minimum length.
    
    :param sequence: A string representing the DNA sequence.
    :param min_length: Minimum length for an ORF (default is 20 nucleotides).
    :return: A list of tuples with (start position, end position, ORF sequence).
    """
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    orfs = []
    seq_len = len(sequence)
    
    for i in range(seq_len - 2):
        codon = sequence[i:i+3]
        if codon == start_codon:  # Found a start codon
            for j in range(i + 3, seq_len - 2, 3):
                stop_codon = sequence[j:j+3]
                if stop_codon in stop_codons:  # Found a stop codon
                    orf = sequence[i:j+3]  # Include the stop codon
                    if len(orf) >= min_length:
                        orfs.append((i, j + 3, orf))
                    break  # Break once the first stop codon is found
    
    return orfs

def write_orfs_to_file(orfs, output_file):
    """
    Writes the ORFs to the output file in a readable format.
    
    :param orfs: List of ORFs as tuples (start, end, orf_sequence).
    :param output_file: Path to the output file.
    """
    with open(output_file, 'w') as out_file:
        for i, (start, end, orf) in enumerate(orfs, start=1):
            out_file.write(f"ORF {i}: Start: {start}, End: {end}, Length: {len(orf)}\n")
            out_file.write(f"{orf}\n\n")

# Example usage:
fasta_file = 'example.fasta'
output_file = 'orfs_output.txt'
min_orf_length = 20  # Minimum ORF length to consider

# Reading the FASTA file
sequences = read_fasta(fasta_file)

# Finding ORFs in each sequence
for header, sequence in sequences.items():
    print(f"Processing {header}...")
    orfs = find_orfs(sequence, min_length=min_orf_length)
    
    if orfs:
        write_orfs_to_file(orfs, output_file)
        print(f"ORFs found and written to {output_file}")
    else:
        print("No ORFs found.")
