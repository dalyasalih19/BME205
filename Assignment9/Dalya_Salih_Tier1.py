import numpy as np
import sys

def read_vcf(vcf_file):
    """Read a VCF file and extract relevant fields."""
    samples = []
    positions = []
    genotypes = {}

    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    header = line.strip().split('\t')
                    samples = header[9:]  # Extract sample names
                    for sample in samples:
                        genotypes[sample] = []
                elif not line.startswith("#"):
                    fields = line.strip().split('\t')
                    pos = int(fields[1])  # Extract genomic position
                    genotype_fields = fields[9:]
                    positions.append(pos)
                    for sample, genotype in zip(samples, genotype_fields):
                        gt = genotype.split(':')[0]  # Extract genotype (e.g., 0|0, 1|1, 0|1)
                        genotypes[sample].append(gt)
    except FileNotFoundError:
        print(f"Error: File '{vcf_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading VCF file: {e}", file=sys.stderr)
        sys.exit(1)

    return positions, genotypes


def viterbi(positions, genotypes, p_ref, error_rate=1 / 1000):
    """Perform Viterbi decoding to identify inbred regions."""
    q_alt = 1 - p_ref
    emission_inbred = {"0|0": 1 - error_rate, "1|1": 1 - error_rate, "0|1": error_rate}
    emission_outbred = {"0|0": 1 - 2 * p_ref * q_alt, "1|1": 1 - 2 * p_ref * q_alt, "0|1": 2 * p_ref * q_alt}

    # Fixed transition probabilities
    transition_inbred_to_outbred = 1 / (1.5 * 10**6)
    transition_outbred_to_inbred = 1 / (4 * 10**6)
    transition_inbred_to_inbred = 1 - transition_inbred_to_outbred
    transition_outbred_to_outbred = 1 - transition_outbred_to_inbred

    results = {}

    for sample, genotype_list in genotypes.items():
        n = len(positions)
        dp = np.full((2, n), -np.inf)  # 0: inbred, 1: outbred
        traceback = np.zeros((2, n), dtype=int)

        # Initialization
        if genotype_list[0] in emission_inbred and genotype_list[0] in emission_outbred:
            dp[0, 0] = np.log(emission_inbred[genotype_list[0]])
            dp[1, 0] = np.log(emission_outbred[genotype_list[0]])

        # Viterbi algorithm
        for i in range(1, n):
            if genotype_list[i] not in emission_inbred or genotype_list[i] not in emission_outbred:
                continue  # Skip invalid genotypes

            for state in range(2):
                if state == 0:  # Inbred
                    trans_probs = [
                        dp[0, i - 1] + np.log(transition_inbred_to_inbred),
                        dp[1, i - 1] + np.log(transition_outbred_to_inbred),
                    ]
                    emission_prob = np.log(emission_inbred[genotype_list[i]])
                else:  # Outbred
                    trans_probs = [
                        dp[0, i - 1] + np.log(transition_inbred_to_outbred),
                        dp[1, i - 1] + np.log(transition_outbred_to_outbred),
                    ]
                    emission_prob = np.log(emission_outbred[genotype_list[i]])

                dp[state, i] = max(trans_probs) + emission_prob
                traceback[state, i] = np.argmax(trans_probs)

        # Traceback to find inbred regions
        current_state = np.argmax(dp[:, -1])
        inbred_regions = []
        end = None

        for i in range(n - 1, -1, -1):
            if current_state == 0:  # Inbred
                if end is None:
                    end = positions[i]
                start = positions[i]
                if i == 0 or traceback[current_state, i] != 0:
                    inbred_regions.append((start, end))
                    end = None
            current_state = traceback[current_state, i]

        results[sample] = sorted(inbred_regions)

    return results


def main():
    if len(sys.argv) != 2:
        print("Usage: python FirstName_LastName_Tier1.py input.vcf", file=sys.stderr)
        return

    vcf_file = sys.argv[1]
    positions, genotypes = read_vcf(vcf_file)

    # Use a fixed reference allele frequency (e.g., 0.5)
    p_ref = 0.5
    results = viterbi(positions, genotypes, p_ref)

    # Print the output in the required tab-delimited format
    print("individual\tstart_position\tstop_position")
    for sample in sorted(results.keys()):  # Sort individuals alphabetically
        regions = sorted(results[sample])  # Sort regions by start position
        for start, stop in regions:
            print(f"{sample}\t{start}\t{stop}")


if __name__ == "__main__":
    main()
