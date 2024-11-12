import argparse
import numpy as np

def simulate_wright_fisher(allele_freq, pop_size, fitness):
    generation_count = 0
    while 0 < allele_freq < 1:
        # Calculate the probability of passing on the allele, adjusted by fitness
        probability = allele_freq * fitness / (allele_freq * fitness + (1 - allele_freq))
        
        # Simulate the next generation allele frequency
        allele_count = np.random.binomial(pop_size, probability)
        allele_freq = allele_count / pop_size
        generation_count += 1
        
    return generation_count, allele_freq

def main():
    parser = argparse.ArgumentParser(description="Simulate allele fixation in a Wright-Fisher model.")
    parser.add_argument('--allele_freq', type=float, required=True, help="Initial allele frequency (0 < freq < 1)")
    parser.add_argument('--pop_size', type=int, required=True, help="Population size (number of individuals)")
    parser.add_argument('--fitness', type=float, required=True, help="Relative fitness of the allele (1=neutral)")
    parser.add_argument('--replicates', type=int, required=True, help="Number of Monte Carlo simulation replicates")
    args = parser.parse_args()
    
    # Variables to track fixation and loss statistics
    fixation_generations = []
    loss_generations = []
    
    for _ in range(args.replicates):
        generation_count, final_freq = simulate_wright_fisher(args.allele_freq, args.pop_size, args.fitness)
        
        # Check if the allele fixed or was lost
        if final_freq == 1:
            fixation_generations.append(generation_count)
        else:
            loss_generations.append(generation_count)
    
    # Calculate statistics for fixation
    if fixation_generations:
        mean_fixation = np.mean(fixation_generations)
        var_fixation = np.var(fixation_generations)
        print(f"Allele was fixed in {mean_fixation:.2f}. Variance: {var_fixation:.2f}")
    else:
        print("Allele was never fixed in any replicate.")
    
    # Calculate statistics for loss
    if loss_generations:
        mean_loss = np.mean(loss_generations)
        var_loss = np.var(loss_generations)
        print(f"Allele was lost in {mean_loss:.2f}. Variance: {var_loss:.2f}")
    else:
        print("Allele was never lost in any replicate.")

if __name__ == "__main__":
    main()
