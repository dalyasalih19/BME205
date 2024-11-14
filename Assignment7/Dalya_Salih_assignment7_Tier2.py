import argparse
import numpy as np
import csv

def read_population_sizes(pop_size_file):
    """
    Reads the population size changes from a TSV file.
    
    Parameters:
    - pop_size_file (str): Path to the TSV file containing generation and population size.
    
    Returns:
    - pop_sizes (list of tuples): Each tuple contains (generation, pop_size).
    """
    pop_sizes = []
    with open(pop_size_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip header
        for row in reader:
            generation, pop_size = int(row[0]), int(row[1])
            pop_sizes.append((generation, pop_size))
    return pop_sizes

def simulate_wright_fisher(allele_freq, fitness, pop_sizes):
    """
    Simulates allele frequency dynamics in a haploid Wright-Fisher population with selection.
    
    Parameters:
    - allele_freq (float): Initial frequency of the allele.
    - fitness (float): Relative fitness of the allele.
    - pop_sizes (list of tuples): List of (generation, population size) tuples.
    
    Returns:
    - generation_count (int): Number of generations until fixation or loss.
    - allele_freq (float): Final allele frequency (either 0 for loss or 1 for fixation).
    """
    generation_count = 0
    current_pop_size = pop_sizes[0][1]
    
    # Loop through generations until fixation or loss
    for gen_idx, (gen_limit, next_pop_size) in enumerate(pop_sizes[1:], start=1):
        while generation_count < gen_limit:
            # Calculate probability of allele transmission based on fitness
            probability = allele_freq * fitness / (allele_freq * fitness + (1 - allele_freq))
            
            # Simulate the next generation's allele count using binomial distribution
            allele_count = np.random.binomial(current_pop_size, probability)
            allele_freq = allele_count / current_pop_size
            generation_count += 1
            
            # Stop if allele reaches fixation or is lost
            if allele_freq == 0 or allele_freq == 1:
                return generation_count, allele_freq
        
        # Update population size for the next interval
        current_pop_size = next_pop_size

    # Continue with final population size if allele hasn't fixed or been lost
    while 0 < allele_freq < 1:
        probability = allele_freq * fitness / (allele_freq * fitness + (1 - allele_freq))
        allele_count = np.random.binomial(current_pop_size, probability)
        allele_freq = allele_count / current_pop_size
        generation_count += 1

    return generation_count, allele_freq

def main():
    parser = argparse.ArgumentParser(description="Simulate allele fixation in a Wright-Fisher model with changing population sizes.")
    parser.add_argument('--allele_freq', type=float, required=True, help="Initial allele frequency (0 < freq < 1)")
    parser.add_argument('--pop_size_file', type=str, required=True, help="TSV file with generation and population size information")
    parser.add_argument('--fitness', type=float, required=True, help="Relative fitness of the allele (1=neutral)")
    parser.add_argument('--replicates', type=int, required=True, help="Number of Monte Carlo simulation replicates")
    args = parser.parse_args()
    
    # Read the population size changes from the file
    pop_sizes = read_population_sizes(args.pop_size_file)
    
    # Lists to collect generation counts for fixation and loss outcomes across replicates
    fixation_generations = []
    loss_generations = []
    
    # Run simulation for the specified number of replicates
    for _ in range(args.replicates):
        generation_count, final_freq = simulate_wright_fisher(args.allele_freq, args.fitness, pop_sizes)
        
        # Append the generation count to either fixation or loss list based on final allele frequency
        if final_freq == 1:
            fixation_generations.append(generation_count)
        else:
            loss_generations.append(generation_count)
    
    # Calculate and output the mean and variance for fixation times if any fixations occurred
    if fixation_generations:
        mean_fixation = np.mean(fixation_generations)
        var_fixation = np.var(fixation_generations)
        print(f"Allele was fixed in {mean_fixation:.2f}. Variance: {var_fixation:.2f}")
    else:
        print("Allele was never fixed in any replicate.")
    
    # Calculate and output the mean and variance for loss times if any losses occurred
    if loss_generations:
        mean_loss = np.mean(loss_generations)
        var_loss = np.var(loss_generations)
        print(f"Allele was lost in {mean_loss:.2f}. Variance: {var_loss:.2f}")
    else:
        print("Allele was never lost in any replicate.")

if __name__ == "__main__":
    main()

