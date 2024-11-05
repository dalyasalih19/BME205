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
    # Define the golden ratio and its conjugate
    phi = (1 + np.sqrt(5)) / 2
    psi = (1 - np.sqrt(5)) / 2

    # Use Binet's formula to get the standard Fibonacci term, F(n)
    F_n = (phi**(n-1) - psi**(n-1)) / np.sqrt(5)

    # Adjust to match the custom starting values by scaling
    if n == 1:
        return a
    elif n == 2:
        return b
    else:
        # Scaling the standard term to fit the custom sequence
        return round(F_n * (b / a))

if __name__ == "__main__":
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    n = int(sys.argv[3])

    print(fibonacci_nth_term(a, b, n))
