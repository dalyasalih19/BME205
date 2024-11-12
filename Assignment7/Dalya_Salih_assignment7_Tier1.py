import argparse
import numpy as np

def simulate_wright_fisher(allele_freq, pop_size, fitness):
    """
    Simulates allele frequency dynamics in a haploid Wright-Fisher population with selection.
    
    Parameters:
    - allele_freq (float): Initial frequency of the allele.
    - pop_size (int): Population size, representing the number of haploid individuals.
    - fitness (float): Relative fitness of the allele. Values >1 indicate positive selection, <1 indicate negative selection.
    
    Returns:
    - generation_count (int): Number of generations until fixation or loss.
    - allele_freq (float): Final allele frequency (either 0 for loss or 1 for fixation).
    """
    generation_count = 0
    
    # Run simulation until the allele either fixes (frequency of 1) or is lost (frequency of 0)
    while 0 < allele_freq < 1:
        # Calculate the probability of passing on the allele, adjusted by fitness
        probability = allele_freq * fitness / (allele_freq * fitness + (1 - allele_freq))
        
        # Simulate the next generation's allele count using binomial distribution
        allele_count = np.random.binomial(pop_size, probability)
        
        # Update allele frequency based on simulated allele count
        allele_freq = allele_count / pop_size
        
        # Increment generation counter
        generation_count += 1
        
    return generation_count, allele_freq

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Simulate allele fixation in a Wright-Fisher model.")
    
    # Add arguments for initial allele frequency, population size, fitness, and number of replicates
    parser.add_argument('--allele_freq', type=float, required=True, help="Initial allele frequency (0 < freq < 1)")
    parser.add_argument('--pop_size', type=int, required=True, help="Population size (number of individuals)")
    parser.add_argument('--fitness', type=float, required=True, help="Relative fitness of the allele (1=neutral)")
    parser.add_argument('--replicates', type=int, required=True, help="Number of Monte Carlo simulation replicates")
    args = parser.parse_args()
    
    # Lists to collect generation counts for fixation and loss outcomes across replicates
    fixation_generations = []
    loss_generations = []
    
    # Run simulation for the specified number of replicates
    for _ in range(args.replicates):
        generation_count, final_freq = simulate_wright_fisher(args.allele_freq, args.pop_size, args.fitness)
        
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

# Run the main function if this script is executed directly
if __name__ == "__main__":
    main()

