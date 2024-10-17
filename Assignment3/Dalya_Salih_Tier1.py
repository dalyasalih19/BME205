import sys
import numpy as np

# Load BED file
def load_bed_file(file_path):
    ranges = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            ranges.append((chromosome, start, end))
    return ranges

# Load genome index
def load_genome_index(file_path):
    chromosome_lengths = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chromosome = parts[0]
            length = int(parts[1])
            chromosome_lengths[chromosome] = length
    return chromosome_lengths

# Merge ranges (Updated)
def merge_ranges_v2(ranges):
    """
    Improved version of merging ranges to ensure correct union of overlapping intervals.
    """
    if not ranges:
        return []
    
    merged = []
    ranges.sort(key=lambda x: (x[0], x[1]))  # Sort by chromosome and start position
    current_chr, current_start, current_end = ranges[0]
    
    for i in range(1, len(ranges)):
        chr_, start, end = ranges[i]
        if chr_ == current_chr and start <= current_end:  # Overlapping or adjacent ranges
            current_end = max(current_end, end)  # Extend the current range
        else:
            merged.append((current_chr, current_start, current_end))  # Add merged range
            current_chr, current_start, current_end = chr_, start, end
    
    merged.append((current_chr, current_start, current_end))  # Append the last range
    return merged

# Calculate overlap using a sweep-line algorithm (Updated)
def calculate_overlap_sweep_v2(merged_a, merged_b):
    total_overlap = 0
    
    i, j = 0, 0
    while i < len(merged_a) and j < len(merged_b):
        chr_a, start_a, end_a = merged_a[i]
        chr_b, start_b, end_b = merged_b[j]
        
        if chr_a == chr_b:
            # Overlap calculation logic
            if end_a < start_b:  # No overlap, move to the next range in Set A
                i += 1
            elif end_b < start_a:  # No overlap, move to the next range in Set B
                j += 1
            else:
                overlap_start = max(start_a, start_b)
                overlap_end = min(end_a, end_b)
                total_overlap += max(0, overlap_end - overlap_start)
                
                if end_a < end_b:
                    i += 1
                else:
                    j += 1
        elif chr_a < chr_b:
            i += 1
        else:
            j += 1
    
    return total_overlap

# Shift ranges for permutation
def shift_ranges(ranges, genome_index):
    shifted = []
    for chr_, start, end in ranges:
        chr_length = genome_index[chr_]
        range_size = end - start
        max_shift = chr_length - range_size
        shift = np.random.randint(0, max_shift)
        shifted.append((chr_, shift, shift + range_size))
    return merge_ranges_v2(shifted)

# Permutation test using the sweep-line algorithm
def permutation_test(merged_a, merged_b, genome_index, n_permutations=10000):
    permuted_overlaps = []
    
    for _ in range(n_permutations):
        permuted_b = shift_ranges(merged_b, genome_index)
        perm_overlap = calculate_overlap_sweep_v2(merged_a, permuted_b)
        permuted_overlaps.append(perm_overlap)
    
    return permuted_overlaps

def main():
    if len(sys.argv) < 4:
        print("Usage: python Dalya_Lastname_Tier1.py <path_to_SetA.bed> <path_to_SetB.bed> <path_to_genome.fa.fai> [n_permutations]")
        sys.exit(1)
    
    set_a_path = sys.argv[1]
    set_b_path = sys.argv[2]
    genome_fai_path = sys.argv[3]
    
    n_permutations = int(sys.argv[4]) if len(sys.argv) > 4 else 10000
    
    # Load files
    set_a = load_bed_file(set_a_path)
    set_b = load_bed_file(set_b_path)
    genome_index = load_genome_index(genome_fai_path)
    
    # Merge ranges (updated version)
    merged_a = merge_ranges_v2(set_a)
    merged_b = merge_ranges_v2(set_b)
    
    # Calculate the observed overlap
    observed_overlap = calculate_overlap_sweep_v2(merged_a, merged_b)
    
    # Run permutation test
    permuted_overlaps = permutation_test(merged_a, merged_b, genome_index, n_permutations)
    
    # Calculate p-value
    p_value = np.sum(np.array(permuted_overlaps) >= observed_overlap) / n_permutations
    
    # Print the result in the required format
    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

if __name__ == "__main__":
    main()
