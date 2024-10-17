import sys
import numpy as np

# Load BED file
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

# Merge ranges (Updated to work per chromosome)
def merge_ranges_v2(ranges):
    merged_ranges = {}
    
    for chromosome, chr_ranges in ranges.items():
        if not chr_ranges:
            continue
        
        merged = []
        chr_ranges.sort(key=lambda x: x[0])  # Sort by start position
        current_start, current_end = chr_ranges[0]
        
        for i in range(1, len(chr_ranges)):
            start, end = chr_ranges[i]
            if start <= current_end:  # Overlapping or adjacent ranges
                current_end = max(current_end, end)  # Extend the current range
            else:
                merged.append((current_start, current_end))  # Add merged range
                current_start, current_end = start, end
        
        merged.append((current_start, current_end))  # Append the last range
        merged_ranges[chromosome] = merged
    
    return merged_ranges

# Calculate overlap using a sweep-line algorithm (Updated for multiple chromosomes)
def calculate_overlap_sweep_v2(merged_a, merged_b):
    total_overlap = 0
    
    for chromosome in merged_a.keys() & merged_b.keys():
        ranges_a = merged_a[chromosome]
        ranges_b = merged_b[chromosome]
        
        i, j = 0, 0
        while i < len(ranges_a) and j < len(ranges_b):
            start_a, end_a = ranges_a[i]
            start_b, end_b = ranges_b[j]
            
            if end_a < start_b:  # No overlap, move to next in Set A
                i += 1
            elif end_b < start_a:  # No overlap, move to next in Set B
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

# Shift ranges for permutation (Chromosome-specific)
def shift_ranges(ranges, genome_index):
    shifted_ranges = {}
    
    for chromosome, chr_ranges in ranges.items():
        chr_length = genome_index[chromosome]
        shifted = []
        
        for start, end in chr_ranges:
            range_size = end - start
            max_shift = chr_length - range_size
            shift = np.random.randint(0, max_shift)
            shifted.append((shift, shift + range_size))
        
        shifted_ranges[chromosome] = merge_ranges_v2({chromosome: shifted})[chromosome]  # Merge shifted ranges
    
    return shifted_ranges

# Permutation test using the sweep-line algorithm (Account for chromosomes)
def permutation_test(merged_a, merged_b, genome_index, n_permutations=10000):
    permuted_overlaps = []
    
    for _ in range(n_permutations):
        permuted_b = shift_ranges(merged_b, genome_index)
        perm_overlap = calculate_overlap_sweep_v2(merged_a, permuted_b)
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
    
    # Merge ranges (updated version)
    merged_a = merge_ranges_v2(set_a)
    merged_b = merge_ranges_v2(set_b)
    
    # Debug: Print the merged ranges for inspection
    print("Merged Set A (sample):", list(merged_a.items())[:2])  # First 2 chromosomes for inspection
    print("Merged Set B (sample):", list(merged_b.items())[:2])  # First 2 chromosomes for inspection
    
    # Calculate the observed overlap
    observed_overlap = calculate_overlap_sweep_v2(merged_a, merged_b)
    
    # Debug: Print the observed overlap
    print(f"Observed overlap: {observed_overlap} bases")
    
    # Run permutation test
    permuted_overlaps = permutation_test(merged_a, merged_b, genome_index, n_permutations)
    
    # Debug: Print a sample of permuted overlaps
    print(f"Sample of permuted overlaps: {permuted_overlaps[:10]}")
    
    # Calculate p-value
    p_value = np.sum(np.array(permuted_overlaps) >= observed_overlap) / n_permutations
    
    # Print the result in the required format
    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

if __name__ == "__main__":
    main()
