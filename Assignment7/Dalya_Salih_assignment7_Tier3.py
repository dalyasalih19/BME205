import argparse
import numpy as np

def simulate_coalescent(pop_size, sample_size, coalescent_event_target):
    """
    Simulates the coalescent process to determine the time until the specified coalescent event.
    
    Parameters:
    - pop_size (int): The population size.
    - sample_size (int): Initial number of lineages (sample size).
    - coalescent_event_target (int): The coalescent event number we are interested in (e.g., 8th event).
    
    Returns:
    - generation_count (int): Number of generations until the specified coalescent event.
    """
    lineages = sample_size
    generation_count = 0
    coalescent_events = 0
    
    while coalescent_events < coalescent_event_target:
        # Calculate the probability of a coalescent event
        prob_coalesce = lineages * (lineages - 1) / (2 * pop_size)
        
        # Calculate the time (in generations) to the next coalescent event
        time_to_event = np.random.geometric(prob_coalesce)
        generation_count += time_to_event
        
        # Perform the coalescent event
        lineages -= 1
        coalescent_events += 1
    
    return generation_count

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Simulate coalescent time to the 8th event in a Wright-Fisher population.")
    parser.add_argument('--pop_size', type=int, required=True, help="Population size (larger background population)")
    parser.add_argument('--sample_size', type=int, required=True, help="Sample size (number of initial lineages)")
    parser.add_argument('--replicates', type=int, required=True, help="Number of Monte Carlo simulation replicates")
    args = parser.parse_args()
    
    # Collect coalescent times for the 8th event across replicates
    coalescent_times = []
    coalescent_event_target = 8
    
    for _ in range(args.replicates):
        generation_count = simulate_coalescent(args.pop_size, args.sample_size, coalescent_event_target)
        coalescent_times.append(generation_count)
    
    # Calculate mean and variance of coalescent times
    mean_coalescent_time = np.mean(coalescent_times)
    var_coalescent_time = np.var(coalescent_times)
    
    # Output results
    print(f"Time to eighth coalescent event: {mean_coalescent_time:.2f}. Variance: {var_coalescent_time:.2f}")

# Run the main function if this script is executed directly
if __name__ == "__main__":
    main()
