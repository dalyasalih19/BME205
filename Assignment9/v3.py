# Filename: FirstName_LastName_Tier1.py
import sys
import numpy as np

def parse_vcf(file_path):
    """Parses a VCF file and extracts relevant information."""
    individuals = []
    positions = []
    genotypes = []
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#CHROM'):
                individuals = line.strip().split('\t')[9:]  # Extract individual names
            elif not line.startswith('#'):
                cols = line.strip().split('\t')
                positions.append(int(cols[1]))  # Position
                genotypes.append(cols[9:])  # Genotypes for individuals
                
    return individuals, positions, np.array(genotypes)

def calculate_emission_probs(genotype, p, e):
    """Calculate emission probabilities for inbred and outbred states."""
    q = 1 - p
    if genotype == '0/0' or genotype == '1/1':  # Homozygous
        inbred_prob = 1 - e
        outbred_prob = 1 - 2 * p * q
    elif genotype == '0/1':  # Heterozygous
        inbred_prob = e
        outbred_prob = 2 * p * q
    else:
        raise ValueError("Unexpected genotype format: {}".format(genotype))
    
    return inbred_prob, outbred_prob

def viterbi_algorithm(positions, genotypes, transition_probs, emission_func, p, e):
    """Runs the Viterbi algorithm to identify the most likely sequence of states."""
    n_positions = len(positions)
    n_states = 2  # Inbred (0) and Outbred (1)
    
    # Transition probabilities
    trans_matrix = np.array([
        [1 - transition_probs['inbred_to_outbred'], transition_probs['inbred_to_outbred']],
        [transition_probs['outbred_to_inbred'], 1 - transition_probs['outbred_to_inbred']]
    ])
    
    log_trans_matrix = np.log(trans_matrix)  # Use log probabilities for numerical stability
    
    # Initialize DP tables
    log_probs = np.full((n_positions, n_states), -np.inf)
    backpointers = np.zeros((n_positions, n_states), dtype=int)
    
    # Initialization
    for state in range(n_states):
        emission_probs = emission_func(genotypes[0], p, e)
        log_probs[0, state] = np.log(emission_probs[state])
    
    # Recursion
    for i in range(1, n_positions):
        for state in range(n_states):
            max_log_prob = -np.inf
            max_state = 0
            for prev_state in range(n_states):
                prob = log_probs[i - 1, prev_state] + log_trans_matrix[prev_state, state]
                if prob > max_log_prob:
                    max_log_prob = prob
                    max_state = prev_state
            log_probs[i, state] = max_log_prob + np.log(emission_func(genotypes[i], p, e)[state])
            backpointers[i, state] = max_state
    
    # Traceback
    states = np.zeros(n_positions, dtype=int)
    states[-1] = np.argmax(log_probs[-1])
    for i in range(n_positions - 2, -1, -1):
        states[i] = backpointers[i + 1, states[i + 1]]
    
    return states

def find_inbred_regions(positions, states, individual):
    """Identifies inbred regions from the Viterbi states."""
    inbred_regions = []
    start = None
    
    for i, state in enumerate(states):
        if state == 0:  # Inbred state
            if start is None:
                start = positions[i]
        elif state == 1 and start is not None:
            inbred_regions.append((individual, start, positions[i - 1]))
            start = None
    
    if start is not None:
        inbred_regions.append((individual, start, positions[-1]))
    
    return inbred_regions

def main(input_file):
    individuals, positions, genotypes = parse_vcf(input_file)
    transition_probs = {
        'inbred_to_outbred': 1 / (1.5 * 10**6),
        'outbred_to_inbred': 1 / (4 * 10**6)
    }
    e = 1 / 1000  # Sequencing error rate
    
    results = []
    for idx, individual in enumerate(individuals):
        states = viterbi_algorithm(
            positions, genotypes[:, idx], transition_probs,
            calculate_emission_probs, p=0.5, e=e
        )
        regions = find_inbred_regions(positions, states, individual)
        results.extend(regions)
    
    # Output results
    print("individual\tstart_position\tstop_position")
    for result in sorted(results, key=lambda x: (x[0], x[1])):
        print(f"{result[0]}\t{result[1]}\t{result[2]}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python FirstName_LastName_Tier1.py input.vcf")
        sys.exit(1)
    main(sys.argv[1])
