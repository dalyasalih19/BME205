import sys
import numpy as np

# Load BED file: This function reads a BED file containing genomic ranges and 
# parses it into a list of tuples where each tuple consists of (chromosome, start, end).
def load_bed_file(file_path):
    ranges = []
    with open(file_path, 'r') as f:
        for line in f:
            # Split each line into its components: chromosome, start, and end.
            parts = line.strip().split()
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            # Append the range as a tuple to the 'ranges' list
            ranges.append((chromosome, start, end))
    return ranges

# Load genome index: This function reads a .fai file (FASTA index file) 
# and stores the length of each chromosome in a dictionary.
def load_genome_index(file_path):
    chromosome_lengths = {}
    with open(file_path, 'r') as f:
        for line in f:
            # Each line contains the chromosome and its length
            parts = line.strip().split()
            chromosome = parts[0]
            length = int(parts[1])
            chromosome_lengths[chromosome] = length
    return chromosome_lengths

# Merge ranges: This function takes a list of genomic ranges and merges any overlapping ranges.
# It sorts the ranges by chromosome and start position and then iterates through them to combine
# overlapping or adjacent intervals.
def merge_ranges_v2(ranges):
    """
    Improved version of merging ranges to ensure correct union of overlapping intervals.
    This function will return the 'union' of the ranges.
    """
    if not ranges:
        return []
    
    merged = []
    # Sort ranges by chromosome and start position
    ranges.sort(key=lambda x: (x[0], x[1]))  
    current_chr, current_start, current_end = ranges[0]
    
    for i in range(1, len(ranges)):
        chr_, start, end = ranges[i]
        # If the next range overlaps or is adjacent, merge it with the current one
        if chr_ == current_chr and start <= current_end:  
            current_end = max(current_end, end)  # Extend the current range
        else:
            # Otherwise, add the current range to the merged list and move to the next range
            merged.append((current_chr, current_start, current_end))  
            current_chr, current_start, current_end = chr_, start, end
    
    # Append the last range after the loop ends
    merged.append((current_chr, current_start, current_end))  
    return merged

# Calculate overlap: This function uses a sweep-line algorithm to calculate the total number 
# of overlapping bases between two merged sets of ranges.
def calculate_overlap_sweep_v2(merged_a, merged_b):
    total_overlap = 0
    i, j = 0, 0  # Initialize indices for both merged range sets
    
    # While there are still ranges to compare
    while i < len(merged_a) and j < len(merged_b):
        chr_a, start_a, end_a = merged_a[i]
        chr_b, start_b, end_b = merged_b[j]
        
        if chr_a == chr_b:
            # If the chromosomes are the same, check for overlap
            if end_a < start_b:  # No overlap, move to the next range in Set A
                i += 1
            elif end_b < start_a:  # No overlap, move to the next range in Set B
                j += 1
            else:
                # Calculate the overlap between the two ranges
                overlap_start = max(start_a, start_b)
                overlap_end = min(end_a, end_b)
                total_overlap += max(0, overlap_end - overlap_start)  # Add to total overlap
                
                # Move to the next range in the set with the earlier end
                if end_a < end_b:
                    i += 1
                else:
                    j += 1
        elif chr_a < chr_b:  # If chromosome in Set A is smaller, move to the next range in Set A
            i += 1
        else:  # Otherwise, move to the next range in Set B
            j += 1
    
    return total_overlap  # Return the total number of overlapping bases

# Shift ranges for permutation: This function randomly shifts ranges within their chromosomes,
# respecting chromosome boundaries, to simulate random sets for the permutation test.
def shift_ranges(ranges, genome_index):
    shifted = []
    for chr_, start, end in ranges:
        chr_length = genome_index[chr_]  # Get the length of the chromosome
        range_size = end - start  # Calculate the size of the range
        max_shift = chr_length - range_size  # Calculate the maximum shift within chromosome
        shift = np.random.randint(0, max_shift)  # Generate a random shift
        # Shift the range and add it to the list
        shifted.append((chr_, shift, shift + range_size))
    # Merge the shifted ranges to ensure no overlapping occurs after shifting
    return merge_ranges_v2(shifted)

# Permutation test: This function performs the permutation test by shifting the ranges in Set B,
# recalculating the overlap, and then comparing the observed overlap with these random overlaps.
def permutation_test(merged_a, merged_b, genome_index, n_permutations=10000):
    permuted_overlaps = []
    
    for _ in range(n_permutations):
        # Shift Set B ranges and calculate the overlap with Set A
        permuted_b = shift_ranges(merged_b, genome_index)
        perm_overlap = calculate_overlap_sweep_v2(merged_a, permuted_b)
        permuted_overlaps.append(perm_overlap)  # Store the permuted overlap
    
    return permuted_overlaps  # Return all permuted overlaps for p-value calculation

def main():
    # Check that the correct number of arguments are passed in
    if len(sys.argv) < 4:
        print("Usage: python Dalya_Salih_Tier1.py <path_to_SetA.bed> <path_to_SetB.bed> <path_to_genome.fa.fai> [n_permutations]")
        sys.exit(1)
    
    # Parse command line arguments for file paths and optional number of permutations
    set_a_path = sys.argv[1]
    set_b_path = sys.argv[2]
    genome_fai_path = sys.argv[3]
    n_permutations = int(sys.argv[4]) if len(sys.argv) > 4 else 10000  # Default to 10,000 permutations
    
    # Load Set A, Set B, and genome index
    set_a = load_bed_file(set_a_path)
    set_b = load_bed_file(set_b_path)
    genome_index = load_genome_index(genome_fai_path)
    
    # Merge ranges in Set A and Set B to form unions
    merged_a = merge_ranges_v2(set_a)
    merged_b = merge_ranges_v2(set_b)
    
    # Calculate the observed overlap between the merged ranges
    observed_overlap = calculate_overlap_sweep_v2(merged_a, merged_b)
    
    # Run the permutation test, shifting Set B and recalculating overlap
    permuted_overlaps = permutation_test(merged_a, merged_b, genome_index, n_permutations)
    
    # Calculate the p-value as the proportion of permuted overlaps greater than or equal to the observed overlap
    p_value = np.sum(np.array(permuted_overlaps) >= observed_overlap) / n_permutations
    
    # Output the result in the required format
    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

if __name__ == "__main__":
    main()

