# DNA Sequence Analyzer

A command-line tool in pure Python that analyzes DNA sequences from FASTA files.

## Features
- GC content calculation
- Reverse complement
- DNA to RNA transcription
- Translation to protein using the standard codon table
- Open Reading Frame (ORF) detection across 3 reading frames

## Usage
```
python analyzer.py sample.fasta --min-orf 5
```

## Sample Output
```
Length: 39 bp
GC content: 56.41%
Reverse complement: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Protein: MAIVMGR*KGAR*
ORF frame 1, 1-24: MAIVMGR
```

## Tools
Python 3

## Future Work
- Reverse strand ORF detection
- Web interface (React + Node.js)