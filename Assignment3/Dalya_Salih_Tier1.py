import sys
import numpy as np

# Interval Tree Node
class IntervalTreeNode:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.left = None
        self.right = None
        self.max_end = end

# Interval Tree Class
class IntervalTree:
    def __init__(self, intervals):
        self.root = self.build_tree(intervals)

    def build_tree(self, intervals):
        if not intervals:
            return None
        
        # Sort intervals by start point
        intervals.sort(key=lambda x: x[0])
        mid = len(intervals) // 2
        
        # Create the root node from the middle of the list
        root = IntervalTreeNode(intervals[mid][0], intervals[mid][1])
        
        # Recursively build the left and right subtrees
        root.left = self.build_tree(intervals[:mid])
        root.right = self.build_tree(intervals[mid + 1:])
        
        # Update max_end to the maximum end value of the node and its children
        root.max_end = max(root.end, 
                           root.left.max_end if root.left else float('-inf'), 
                           root.right.max_end if root.right else float('-inf'))
        
        return root

    def query(self, node, start, end):
        if node is None:
            return []
        
        result = []
        
        # Check if the current interval overlaps with the query interval
        if node.start <= end and start <= node.end:
            result.append((node.start, node.end))
        
        # Recurse into the left subtree if it might contain overlapping intervals
        if node.left and node.left.max_end >= start:
            result.extend(self.query(node.left, start, end))
        
        # Recurse into the right subtree
        result.extend(self.query(node.right, start, end))
        
        return result

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

# Merge ranges
def merge_ranges(ranges):
    merged = []
    ranges.sort(key=lambda x: (x[0], x[1]))
    current_chr, current_start, current_end = ranges[0]
    
    for i in range(1, len(ranges)):
        chr_, start, end = ranges[i]
        if chr_ == current_chr and start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_chr, current_start, current_end))
            current_chr, current_start, current_end = chr_, start, end
    
    merged.append((current_chr, current_start, current_end))
    return merged

# Calculate overlap using the interval tree
def calculate_overlap_with_tree(merged_a, merged_b):
    total_overlap = 0
    chr_to_intervals = {}
    
    # Build interval trees for each chromosome
    for chr_, start, end in merged_a:
        if chr_ not in chr_to_intervals:
            chr_to_intervals[chr_] = []
        chr_to_intervals[chr_].append((start, end))
    
    trees = {chr_: IntervalTree(intervals) for chr_, intervals in chr_to_intervals.items()}
    
    # Query the tree for overlaps with Set B
    for chr_, start_b, end_b in merged_b:
        if chr_ in trees:
            overlaps = trees[chr_].query(trees[chr_].root, start_b, end_b)
            for start_a, end_a in overlaps:
                overlap_start = max(start_a, start_b)
                overlap_end = min(end_a, end_b)
                total_overlap += max(0, overlap_end - overlap_start)
    
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
    return merge_ranges(shifted)

# Permutation test using the interval tree
def permutation_test(merged_a, merged_b, genome_index, n_permutations=10000):
    permuted_overlaps = []
    
    for _ in range(n_permutations):
        permuted_b = shift_ranges(merged_b, genome_index)
        perm_overlap = calculate_overlap_with_tree(merged_a, permuted_b)
        permuted_overlaps.append(perm_overlap)
    
    return permuted_overlaps

def main():
    if len(sys.argv) < 4:
        print("Usage: python Dalya_Salih_Tier1.py <path_to_SetA.bed> <path_to_SetB.bed> <path_to_genome.fa.fai> [n_permutations]")
        sys.exit(1)
    
    set_a_path = sys.argv[1]
    set_b_path = sys.argv[2]
    genome_fai_path = sys.argv[3]
    
    n_permutations = int(sys.argv[4]) if len(sys.argv) > 4 else 10000
    
    set_a = load_bed_file(set_a_path)
    set_b = load_bed_file(set_b_path)
    genome_index = load_genome_index(genome_fai_path)
    
    merged_a = merge_ranges(set_a)
    merged_b = merge_ranges(set_b)
    
    observed_overlap = calculate_overlap_with_tree(merged_a, merged_b)
    
    permuted_overlaps = permutation_test(merged_a, merged_b, genome_index, n_permutations)
    
    p_value = np.sum(np.array(permuted_overlaps) >= observed_overlap) / n_permutations
    
    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

if __name__ == "__main__":
    main()

