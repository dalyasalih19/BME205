import numpy as np
import sys

def fibonacci_custom_eigen(a, b, n):
    # Define the Fibonacci transformation matrix, M
    # This matrix transforms a vector representing consecutive Fibonacci terms
    # to the next terms in the sequence.
    M = np.array([[1, 1],
                  [1, 0]])

    # Calculate eigenvalues and eigenvectors of the matrix M
    # The eigenvalues (alpha and beta) of M are the roots of its characteristic polynomial,
    # which relate to the golden ratio. These values help us derive the closed-form formula
    # for the Fibonacci sequence.
    eigenvalues, eigenvectors = np.linalg.eig(M)
    
    # Extract the eigenvalues
    # alpha is the positive root (golden ratio) and beta is the negative root.
    alpha, beta = eigenvalues  # alpha > beta, corresponds to phi and -1/phi

    # The matrix M has two eigenvectors, corresponding to eigenvalues alpha and beta.
    # v_alpha and v_beta are vectors that indicate directions along which the transformation
    # scales by alpha or beta, respectively.
    v_alpha = eigenvectors[:, 0]  # eigenvector associated with alpha (golden ratio)
    v_beta = eigenvectors[:, 1]   # eigenvector associated with beta

    # Calculate coefficients c1 and c2 for the linear combination of eigenvectors
    # that represent the initial state of the sequence.
    # We want the linear combination of v_alpha and v_beta that equals the starting
    # values [F(2), F(1)], where F(2) = b and F(1) = a.
    # This setup yields the equation: initial_vector = c1 * v_alpha + c2 * v_beta.
    initial_vector = np.array([b, a])  # Starting vector with custom values [F(2), F(1)]

    # Form the matrix A using v_alpha and v_beta as columns so that we can solve for
    # c1 and c2 in the equation A * [c1, c2] = initial_vector.
    A = np.column_stack((v_alpha, v_beta))
    c1, c2 = np.linalg.solve(A, initial_vector)

    # Compute the n-th term of the custom sequence using the closed-form solution.
    # The n-th term is derived by raising the eigenvalues (alpha and beta) to the (n-2) power,
    # then scaling by c1 and c2 and the respective eigenvector components.
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


