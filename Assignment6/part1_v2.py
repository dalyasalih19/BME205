import numpy as np
import sys

def fibonacci_custom_eigen(a, b, n):
    # Define the Fibonacci matrix
    M = np.array([[1, 1],
                  [1, 0]])

    # Calculate eigenvalues and eigenvectors of the matrix
    eigenvalues, eigenvectors = np.linalg.eig(M)
    
    # Extract the eigenvalues
    alpha, beta = eigenvalues  # alpha > beta, corresponds to phi and -1/phi

    # Eigenvector corresponding to alpha (golden ratio) and beta
    v_alpha = eigenvectors[:, 0]  # eigenvector associated with alpha
    v_beta = eigenvectors[:, 1]   # eigenvector associated with beta

    # Calculate coefficients c1 and c2 to match initial terms a and b
    # Set up system to solve: initial_vector = c1 * v_alpha + c2 * v_beta
    initial_vector = np.array([b, a])  # [F(2), F(1)] as starting values
    A = np.column_stack((v_alpha, v_beta))
    c1, c2 = np.linalg.solve(A, initial_vector)

    # Compute the n-th term using eigenvalues raised to power (n-2)
    term_n = c1 * (alpha ** (n - 2)) * v_alpha[0] + c2 * (beta ** (n - 2)) * v_beta[0]

    # Return the term as an integer (rounded)
    return int(round(term_n))

if __name__ == "__main__":
    # Parse command-line arguments
    start_a = int(sys.argv[1])
    start_b = int(sys.argv[2])
    fNumber = int(sys.argv[3])

    # Calculate and print the n-th Fibonacci term for custom starting values
    print(fibonacci_custom_eigen(start_a, start_b, fNumber))

