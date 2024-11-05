import numpy as np
import sys

def fibonacci_custom_sequence(a, b, n):
    # Generate the indices for Fibonacci terms
    indices = np.arange(1, n + 1)
    length_indices = len(indices)

    # Constants for Binet's formula
    sqrt_five = np.sqrt(5) 
    alpha = (1 + sqrt_five) / 2
    beta = (1 - sqrt_five) / 2

    # Calculate the standard Fibonacci sequence terms using Binet's formula
    Fn_standard = np.rint(((alpha ** indices) - (beta ** indices)) / sqrt_five)

    # Adjust the terms based on custom starting values
    Fn_custom = np.empty(length_indices)
    Fn_custom[0] = a
    Fn_custom[1] = b
    for i in range(2, length_indices):
        Fn_custom[i] = Fn_standard[i - 1] * (b / Fn_standard[1]) + Fn_standard[i - 2] * (a / Fn_standard[0])

    # Convert to integers for readability
    Fn_custom = Fn_custom.astype(int)
    
    # Output the n-th term
    return Fn_custom

if __name__ == "__main__":
    # Command-line arguments for custom starting values and the desired term count
    start_a = int(sys.argv[1])
    start_b = int(sys.argv[2])
    fNumber = int(sys.argv[3])

    # Calculate the custom sequence
    fibonacci_sequence = fibonacci_custom_sequence(start_a, start_b, fNumber)

    # Output the sequence and the n-th term
    print("The first {} numbers of the Fibonacci series with custom starting values are: {}.".format(fNumber, fibonacci_sequence))
    print({fibonacci_sequence[-1]})


