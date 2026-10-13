import argparse

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
CODON_TABLE = dict(zip(CODONS, AMINO))


def clean_sequence(seq):
    seq = "".join(seq.split()).upper()
    if not set(seq) <= set("ACGT"):
        raise ValueError("Sequence e shudhu A, C, G, T thakte pare")
    return seq


def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def transcribe(seq):
    return seq.replace("T", "U")


def translate(seq):
    protein = ""
    for i in range(0, len(seq) - 2, 3):
        protein += CODON_TABLE[seq[i:i+3]]
    return protein


def find_orfs(seq, min_length=30):
    orfs = []
    for frame in range(3):
        i = frame
        while i + 3 <= len(seq):
            if seq[i:i+3] == "ATG":
                protein = ""
                for j in range(i, len(seq) - 2, 3):
                    aa = CODON_TABLE[seq[j:j+3]]
                    if aa == "*":
                        if len(protein) >= min_length:
                            orfs.append((i + 1, j + 3, frame + 1, protein))
                        i = j
                        break
                    protein += aa
            i += 3
    return orfs


def read_fasta(path):
    with open(path) as f:
        lines = [line.strip() for line in f if not line.startswith(">")]
    return "".join(lines)


def main():
    parser = argparse.ArgumentParser(description="DNA Sequence Analyzer")
    parser.add_argument("fasta", help="FASTA file er path")
    parser.add_argument("--min-orf", type=int, default=30)
    args = parser.parse_args()

    seq = clean_sequence(read_fasta(args.fasta))
    print(f"Length: {len(seq)} bp")
    print(f"GC content: {gc_content(seq):.2f}%")
    print(f"Reverse complement: {reverse_complement(seq)}")
    print(f"RNA: {transcribe(seq)}")
    print(f"Protein: {translate(seq)}")
    for start, end, frame, protein in find_orfs(seq, args.min_orf):
        print(f"ORF frame {frame}, {start}-{end}: {protein}")


if __name__ == "__main__":
    main()