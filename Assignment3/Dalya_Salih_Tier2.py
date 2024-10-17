import sys
import numpy as np

# Load BED file with chromosome-specific ranges
def load_bed_file(file_path):
    ranges = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if chromosome not in ranges:
                ranges[chromosome] = []
            ranges[chromosome].append((start, end))
    return ranges

# Load genome index for chromosome lengths
def load_genome_index(file_path):
    chromosome_lengths = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chromosome = parts[0]
            length = int(parts[1])
            chromosome_lengths[chromosome] = length
    return chromosome_lengths

# Merge ranges within a set for a specific chromosome
def merge_ranges(ranges):
    if not ranges:
        return []
    
    merged = []
    ranges.sort()  # Sort by start position
    current_start, current_end = ranges[0]
    
    for start, end in ranges[1:]:
        if start <= current_end:  # Overlapping or adjacent ranges
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    
    merged.append((current_start, current_end))  # Append the last range
    return merged

# Calculate overlap using a sweep-line algorithm for specific chromosome
def calculate_overlap(merged_a, merged_b):
    total_overlap = 0
    i, j = 0, 0
    
    while i < len(merged_a) and j < len(merged_b):
        start_a, end_a = merged_a[i]
        start_b, end_b = merged_b[j]
        
        if end_a < start_b:
            i += 1
        elif end_b < start_a:
            j += 1
        else:
            overlap_start = max(start_a, start_b)
            overlap_end = min(end_a, end_b)
            total_overlap += max(0, overlap_end - overlap_start)
            
            if end_a < end_b:
                i += 1
            else:
                j += 1
    
    return total_overlap

# Shift ranges within a specific chromosome for permutation
def shift_ranges(ranges, chromosome_length):
    shifted = []
    for start, end in ranges:
        range_size = end - start
        max_shift = chromosome_length - range_size
        shift = np.random.randint(0, max_shift)
        shifted.append((shift, shift + range_size))
    return merge_ranges(shifted)

# Permutation test considering chromosome structure
def permutation_test(merged_a, merged_b, genome_index, n_permutations=10000):
    permuted_overlaps = []
    
    for _ in range(n_permutations):
        permuted_b_all = {}
        for chr_ in merged_b:
            if chr_ in genome_index:
                permuted_b_all[chr_] = shift_ranges(merged_b[chr_], genome_index[chr_])
        
        # Calculate total overlap across all chromosomes
        perm_overlap = 0
        for chr_ in merged_a:
            if chr_ in permuted_b_all:
                perm_overlap += calculate_overlap(merged_a[chr_], permuted_b_all[chr_])
        
        permuted_overlaps.append(perm_overlap)
    
    return permuted_overlaps

def main():
    if len(sys.argv) < 4:
        print("Usage: python Dalya_Lastname_Tier2.py <path_to_SetA.bed> <path_to_SetB.bed> <path_to_genome.fa.fai> [n_permutations]")
        sys.exit(1)
    
    set_a_path = sys.argv[1]
    set_b_path = sys.argv[2]
    genome_fai_path = sys.argv[3]
    
    n_permutations = int(sys.argv[4]) if len(sys.argv) > 4 else 10000
    
    # Load files
    set_a = load_bed_file(set_a_path)
    set_b = load_bed_file(set_b_path)
    genome_index = load_genome_index(genome_fai_path)
    
    # Merge ranges for each chromosome
    merged_a = {chr_: merge_ranges(ranges) for chr_, ranges in set_a.items()}
    merged_b = {chr_: merge_ranges(ranges) for chr_, ranges in set_b.items()}
    
    # Calculate the observed overlap
    observed_overlap = 0
    for chr_ in merged_a:
        if chr_ in merged_b:
            observed_overlap += calculate_overlap(merged_a[chr_], merged_b[chr_])
    
    # Run permutation test
    permuted_overlaps = permutation_test(merged_a, merged_b, genome_index, n_permutations)
    
    # Calculate p-value
    p_value = np.sum(np.array(permuted_overlaps) >= observed_overlap) / n_permutations
    
    # Print the result in the required format
    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

if __name__ == "__main__":
    main()

