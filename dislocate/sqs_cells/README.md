# Special Quasirandom Structures (SQS)

This directory contains Special Quasirandom Structure (SQS) files of hcp, fcc, bcc, and omega phase for binary and ternary alloy systems. All structures contains 40 atoms.

## Directory Structure

```
sqs_cells/
├── binary/              # Binary alloy SQS structures
│   ├── bcc_A0.10_B0.90/ # BCC structure, 10% A, 90% B
│   ├── bcc_A0.20_B0.80/ # BCC structure, 20% A, 80% B
│   ├── fcc_A0.10_B0.90/ # FCC structure, 10% A, 90% B
│   ├── fcc_A0.20_B0.80/ # FCC structure, 20% A, 80% B
│   ├── hcp_A0.10_B0.90/ # HCP structure, 10% A, 90% B
│   ├── hcp_A0.20_B0.80/ # HCP structure, 20% A, 80% B
│   ├── omega_A0.10_B0.90/ # Omega structure, 10% A, 90% B
│   └── omega_A0.20_B0.80/ # Omega structure, 20% A, 80% B
└── ternary/             # Ternary alloy SQS structures
    ├── bcc_A0.125_B0.125_C0.750/ # BCC structure, 12.5% A, 12.5% B, 75% C
    ├── fcc_A0.125_B0.125_C0.750/ # FCC structure, 12.5% A, 12.5% B, 75% C
    ├── hcp_A0.125_B0.125_C0.750/ # HCP structure, 12.5% A, 12.5% B, 75% C
```

## File Formats

### Output Files
Each directory contains multiple output files from SQS generation:

- **`bestsqs1.out`** to **`bestsqs5.out`**: Best SQS structures found
- **`bestcorr1.out`** to **`bestcorr5.out`**: Correlation function data
- **`clusters.out`**: Cluster correlation information

## Crystal Systems

### Binary Alloys
- **BCC**: Body-centered cubic structures
- **FCC**: Face-centered cubic structures  
- **HCP**: Hexagonal close-packed structures
- **Omega**: Omega phase structures

### Ternary Alloys
- **BCC**: Body-centered cubic ternary alloys
- **FCC**: Face-centered cubic ternary alloys
- **HCP**: Hexagonal close-packed ternary alloys
- **Omega**: Omega phase ternary alloys

## Compositions

### Binary Compositions
- **A0.10_B0.90**: 10% A, 90% B
- **A0.20_B0.80**: 20% A, 80% B
- **A0.30_B0.70**: 30% A, 70% B
- **A0.40_B0.60**: 40% A, 60% B
- **A0.50_B0.50**: 50% A, 50% B

### Ternary Compositions
- **A0.125_B0.125_C0.750**: 12.5% A, 12.5% B, 75% C
- **A0.250_B0.250_C0.500**: 25% A, 25% B, 50% C
- **A0.333_B0.333_C0.333**: 33.3% A, 33.3% B, 33.3% C

## References

- **ATAT**: Alloy Theoretic Automated Toolkit
- **SQS**: Special Quasirandom Structures
- **Cluster Expansion**: Statistical mechanics of alloys 