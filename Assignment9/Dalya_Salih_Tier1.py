import numpy as np
import sys

def read_vcf(vcf_file):
    """Read a VCF file and extract relevant fields."""
    samples = []
    positions = []
    genotypes = {}

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().split('\t')
                samples = header[9:]  # Sample names
                for sample in samples:
                    genotypes[sample] = []
            elif not line.startswith("#"):
                fields = line.strip().split('\t')
                pos = int(fields[1])  # Position
                genotype_fields = fields[9:]
                positions.append(pos)
                for sample, genotype in zip(samples, genotype_fields):
                    gt = genotype.split(':')[0]  # Extract genotype
                    genotypes[sample].append(gt)

    return positions, genotypes


def viterbi(positions, genotypes, p_ref, error_rate=1/1000):
    """Perform Viterbi decoding to identify inbred regions."""
    q_alt = 1 - p_ref
    emission_inbred = {"0|0": 1 - error_rate, "1|1": 1 - error_rate, "0|1": error_rate}
    emission_outbred = {"0|0": 1 - 2 * p_ref * q_alt, "1|1": 1 - 2 * p_ref * q_alt, "0|1": 2 * p_ref * q_alt}

    transition_inbred_to_outbred = 1 / (1.5 * 10**6)
    transition_outbred_to_inbred = 1 / (4 * 10**6)
    transition_inbred_to_inbred = 1 - transition_inbred_to_outbred
    transition_outbred_to_outbred = 1 - transition_outbred_to_inbred

    results = {}

    for sample, genotype_list in genotypes.items():
        n = len(positions)
        dp = np.zeros((2, n))  # 0: inbred, 1: outbred
        traceback = np.zeros((2, n), dtype=int)

        # Initialization
        if genotype_list[0] in emission_inbred and genotype_list[0] in emission_outbred:
            dp[0, 0] = np.log(emission_inbred[genotype_list[0]])
            dp[1, 0] = np.log(emission_outbred[genotype_list[0]])

        # Viterbi algorithm
        for i in range(1, n):
            delta_pos = positions[i] - positions[i - 1]

            for state in range(2):
                if genotype_list[i] not in emission_inbred or genotype_list[i] not in emission_outbred:
                    continue

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
        end = positions[-1]

        for i in range(n - 1, -1, -1):
            if current_state == 0:  # Inbred
                start = positions[i]
                if i == 0 or traceback[current_state, i] != 0:
                    inbred_regions.append((start, end))
                    end = positions[i - 1] if i > 0 else start

            current_state = traceback[current_state, i]

        results[sample] = sorted(inbred_regions)

    return results


def main():
    if len(sys.argv) != 2:
        print("Usage: python Dalya_Salih_Tier1.py synthetic_population.vcf")
        return

    vcf_file = sys.argv[1]
    positions, genotypes = read_vcf(vcf_file)

    # Assuming a fixed reference allele frequency (e.g., 0.5)
    p_ref = 0.5
    results = viterbi(positions, genotypes, p_ref)

    print("individual\tstart_position\tstop_position")
    for sample, regions in sorted(results.items()):
        for start, stop in regions:
            print(f"{sample}\t{start}\t{stop}")


if __name__ == "__main__":
    main()
