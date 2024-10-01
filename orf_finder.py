Last login: Tue Oct  1 11:20:11 on ttys000
(base) admin.dsalih@ITS-dsalih-lt ~ % conda activate dalya_env  
(dalya_env) admin.dsalih@ITS-dsalih-lt ~ % vim example.fasta
(dalya_env) admin.dsalih@ITS-dsalih-lt ~ % ls
Applications				Desktop					Library					Pictures				miniconda3
BigFix					Documents				Movies					Public
Creative Cloud Files			Downloads				Music					Salih_Dalya_BME163_Assignment_Week2.py
(dalya_env) admin.dsalih@ITS-dsalih-lt ~ % cd Desktop 
(dalya_env) admin.dsalih@ITS-dsalih-lt Desktop % ls
BME163						Personal					Screenshot 2024-09-06 at 2.02.44 PM.png		Screenshot 2024-09-30 at 12.25.58 PM.png
CSE101						Screenshot 2024-07-23 at 2.19.27 PM.png		Screenshot 2024-09-06 at 2.03.14 PM.png		Screenshot 2024-09-30 at 12.26.02 PM.png
CSE182						Screenshot 2024-08-13 at 4.39.27 PM.png		Screenshot 2024-09-19 at 9.27.17 AM.png		advent_calendaR
Lab						Screenshot 2024-08-26 at 1.45.18 PM.png		Screenshot 2024-09-26 at 3.11.31 PM.png
(dalya_env) admin.dsalih@ITS-dsalih-lt Desktop % ls
BME163						Lab						Screenshot 2024-08-26 at 1.45.18 PM.png		Screenshot 2024-09-26 at 3.11.31 PM.png
BME205						Personal					Screenshot 2024-09-06 at 2.02.44 PM.png		Screenshot 2024-09-30 at 12.25.58 PM.png
CSE101						Screenshot 2024-07-23 at 2.19.27 PM.png		Screenshot 2024-09-06 at 2.03.14 PM.png		Screenshot 2024-09-30 at 12.26.02 PM.png
CSE182						Screenshot 2024-08-13 at 4.39.27 PM.png		Screenshot 2024-09-19 at 9.27.17 AM.png		advent_calendaR
(dalya_env) admin.dsalih@ITS-dsalih-lt Desktop % cd BME205
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % vim example.fasta                                 
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % vim orf_finder.py
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % python orf_finder.py
Processing sequence_1...
ORFs found and written to orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % ls
example.fasta	orf_finder.py	orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % less orfs_output.txt 
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % less orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % less example.fasta  
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % ls
example.fasta	orf_finder.py	orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % rm orf_finder.py 
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % vim orf_finder.py
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % python orf_finder.py
Processing sequence_1...
ORFs found and written to orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % less orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % ls
example.fasta	orf_finder.py	orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % rm example.fasta 
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % vim example.fasta
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % python orf_finder.py
Processing sequence_1...
ORFs found and written to orfs_output.txt
Processing sequence_2...
ORFs found and written to orfs_output.txt
Processing sequence_3...
ORFs found and written to orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % less orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % rm orf_finder.py 
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % vim orf_finder.py
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % python orf_finder.py
Processing sequence_1...
ORFs found and written to orfs_output.txt
Processing sequence_2...
ORFs found and written to orfs_output.txt
Processing sequence_3...
ORFs found and written to orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % less orfs_output.txt
(dalya_env) admin.dsalih@ITS-dsalih-lt BME205 % less orf_finder.py 



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
    specified minimum length, including overlapping ORFs.
    
    :param sequence: A string representing the DNA sequence.
    :param min_length: Minimum length for an ORF in nucleotides (default is 20 nucleotides).
    :return: A list of tuples with (start position, end position, ORF sequence).
    """
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    orfs = []
    seq_len = len(sequence)
    
    # Search for all start codons and then find the first stop codon after each start
    for i in range(seq_len - 2):
        codon = sequence[i:i+3]
        if codon == start_codon:  # Found a start codon
            for j in range(i + 3, seq_len - 2, 3):
                stop_codon = sequence[j:j+3]
                if stop_codon in stop_codons:  # Found a stop codon
                    orf = sequence[i:j+3]  # Include the stop codon
                    if len(orf) >= min_length:
                        orfs.append((i, j + 3, orf))
                    # Do not break here; continue searching for other stop codons in case of overlapping ORFs
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
min_orf_length = 20  # Minimum ORF length to consider (in nucleotides)

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
        
