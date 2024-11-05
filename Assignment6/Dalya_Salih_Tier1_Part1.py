# import numpy as np
# import sys

# def fibonacci_nth_term(a, b, n):
#     # Define the transformation matrix for Fibonacci sequence
#     matrix = np.array([[1, 1], [1, 0]])
    
#     # Find eigenvalues and eigenvectors
#     eigenvalues, eigenvectors = np.linalg.eig(matrix)
    
#     # Use the eigenvalues and eigenvectors to derive the explicit formula for Fibonacci-like sequences
#     phi = (1 + np.sqrt(5)) / 2
#     psi = (1 - np.sqrt(5)) / 2

#     # Calculate the n-th term using the closed-form formula with eigenvalues
#     term_n = (a * (phi**(n - 1)) - b * (psi**(n - 1))) / np.sqrt(5)
#     return round(term_n)

# if __name__ == "__main__":
#     a = int(sys.argv[1])
#     b = int(sys.argv[2])
#     n = int(sys.argv[3])

#     print(fibonacci_nth_term(a, b, n))


import numpy as np
import sys

def fibonacci_nth_term(a, b, n):
    # Define the golden ratio components
    sqrt_five = np.sqrt(5)
    alpha = (1 + sqrt_five) / 2
    beta = (1 - sqrt_five) / 2

    # Use Binet's formula to calculate the n-th Fibonacci term
    # The formula assumes the standard Fibonacci sequence starting at F(1) = 1 and F(2) = 1
    # We calculate F(n-2) + F(n-1) for terms beyond the first two
    if n == 1:
        return a
    elif n == 2:
        return b
    else:
        # Calculate the first two terms based on standard Fibonacci values
        F1 = a
        F2 = b
        # Compute the n-th term using Binet's formula adjusted for starting values
        Fn = np.rint(((alpha ** (n-1)) - (beta ** (n-1))) / sqrt_five).astype(int)
        # Scale Fn to match the custom starting sequence
        return round(F1 + (Fn - 1) * (F2 / F1))

if __name__ == "__main__":
    # Input parameters
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    n = int(sys.argv[3])

    # Output the n-th term
    print(fibonacci_nth_term(a, b, n))

