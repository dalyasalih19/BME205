# Filename: FirstName_LastName_Tier1.py

import sys
import numpy as np

# Constants
INBRED_TO_OUTBRED = 1 / (1.5 * 10**6)
OUTBRED_TO_INBRED = 1 / (4 * 10**6)
ERROR_RATE = 1 / 1000  # Sequencing error rate
DEFAULT_ALLELE_FREQ = 0.5  # Assumed reference allele frequency

def parse_vcf(vcf_file):
    """Parse a VCF file and extract positions, genotypes, and individuals."""
    positions = []
    genotypes = []
    individuals = []
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith("#CHROM"):  # Header line with individual names
                header = line.strip().split('\t')
                individuals = header[9:]  # Extract individual names from the VCF header
            elif not line.startswith("#"):
                fields = line.strip().split('\t')
                positions.append(int(fields[1]))  # Position of the variant
                genotypes.append(fields[9:])  # Genotypes for all individuals at this position

    return np.array(positions), np.array(genotypes), individuals

def calculate_emission_probs(genotype, p):
    """Calculate emission probabilities for inbred and outbred states."""
    q = 1 - p
    genotype = genotype.replace('|', '/')  # Treat phased genotypes as unphased
    if genotype in {'0/0', '1/1'}:  # Homozygous
        inbred_prob = 1 - ERROR_RATE
        outbred_prob = 1 - 2 * p * q
    elif genotype in {'0/1', '1/0'}:  # Heterozygous
        inbred_prob = ERROR_RATE
        outbred_prob = 2 * p * q
    else:
        raise ValueError(f"Unexpected genotype format: {genotype}")
    
    return inbred_prob, outbred_prob

def viterbi(positions, genotypes, transition_probs, p):
    """Perform the Viterbi algorithm to find the most likely sequence of inbred states."""
    n_positions = len(positions)
    n_states = 2  # 0: Inbred, 1: Outbred
    log_probs = np.full((n_positions, n_states), -np.inf)
    backpointers = np.zeros((n_positions, n_states), dtype=int)
    
    # Transition matrix in log space
    trans_matrix = np.array([
        [1 - transition_probs['inbred_to_outbred'], transition_probs['inbred_to_outbred']],
        [transition_probs['outbred_to_inbred'], 1 - transition_probs['outbred_to_inbred']]
    ])
    log_trans_matrix = np.log(trans_matrix)
    
    # Initialization
    for state in range(n_states):
        emission_probs = calculate_emission_probs(genotypes[0], p)
        log_probs[0, state] = np.log(emission_probs[state])
    
    # Forward pass
    for i in range(1, n_positions):
        for state in range(n_states):
            max_prob, max_state = max(
                (log_probs[i-1, prev_state] + log_trans_matrix[prev_state, state], prev_state)
                for prev_state in range(n_states)
            )
            emission_probs = calculate_emission_probs(genotypes[i], p)
            log_probs[i, state] = max_prob + np.log(emission_probs[state])
            backpointers[i, state] = max_state
    
    # Traceback
    states = np.zeros(n_positions, dtype=int)
    states[-1] = np.argmax(log_probs[-1])
    for i in range(n_positions - 2, -1, -1):
        states[i] = backpointers[i + 1, states[i + 1]]
    
    return states

def find_inbred_regions(positions, states):
    """Identify inbred regions from the Viterbi states."""
    regions = []
    start = None
    for i, state in enumerate(states):
        if state == 0:  # Inbred state
            if start is None:
                start = positions[i]
        elif state == 1 and start is not None:
            regions.append((start, positions[i - 1]))
            start = None
    if start is not None:
        regions.append((start, positions[-1]))
    return regions

def sort_key(result):
    """Sort key function for ordering by individual and start position."""
    name, start, stop = result
    numeric_part = int(''.join(filter(str.isdigit, name)) or 0)
    return (numeric_part, start)

def main(vcf_file):
    positions, genotypes, individuals = parse_vcf(vcf_file)
    transition_probs = {
        'inbred_to_outbred': INBRED_TO_OUTBRED,
        'outbred_to_inbred': OUTBRED_TO_INBRED
    }
    p = DEFAULT_ALLELE_FREQ  # Assume equal reference/alternative frequencies for now
    
    results = []
    for i, individual in enumerate(individuals):
        individual_genotypes = genotypes[:, i]
        states = viterbi(positions, individual_genotypes, transition_probs, p)
        regions = find_inbred_regions(positions, states)
        for start, end in regions:
            results.append((individual, start, end))
    
    # Sort results by individual name and start position
    results.sort(key=sort_key)
    
    # Output results in the specified format
    print("individual\tstart_position\tstop_position")
    for individual, start, stop in results:
        print(f"{individual}\t{start}\t{stop}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python FirstName_LastName_Tier1.py input.vcf")
        sys.exit(1)
    main(sys.argv[1])

