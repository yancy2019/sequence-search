# sequence-search
#blast proline rich region
#Usage:
#blast_local_proline_number.py [-h] [--region_length REGION_LENGTH] [--proline_threshold PROLINE_THRESHOLD] input_file output_file

def count_proline_in_region(sequence, region_length=500, proline_threshold=200):
    for start in range(len(sequence) - region_length + 1):
        region = sequence[start:start + region_length].upper()
        if region.count('P') > proline_threshold:
            return True
    return False

def filter_sequences(input_file, output_file, region_length=500, proline_threshold=200):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_header = ""
        current_sequence = ""
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if current_sequence and count_proline_in_region(current_sequence, region_length, proline_threshold):
                    outfile.write(">{header}\n{sequence}\n".format(header=current_header, sequence=current_sequence))
                current_header = line[1:]
                current_sequence = ""
            else:
                current_sequence += line

        if current_sequence and count_proline_in_region(current_sequence, region_length, proline_threshold):
            outfile.write(">{header}\n{sequence}\n".format(header=current_header, sequence=current_sequence))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Filter protein sequences with more than a specified number of 'P' residues in any local region of a given length.")
    parser.add_argument("input_file", help="Path to the input FASTA file.")
    parser.add_argument("output_file", help="Path to the output filtered FASTA file.")
    parser.add_argument("--region_length", type=int, default=500, help="Length of the local region (default: 500).")
    parser.add_argument("--proline_threshold", type=int, default=200, help="Threshold for the number of 'P' residues (default: 200).")
    args = parser.parse_args()

    filter_sequences(args.input_file, args.output_file, args.region_length, args.proline_threshold) 
