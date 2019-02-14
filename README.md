# BBS
Bio Bank Simulator. A simulation software for generating phenotypes from UK biobank data

Biobank Simulation Tool
Usage: BBS [options]
Options:
    --input   | -i    Input file prefix
    --out     | -o    Output file prefix
    --num     | -n    Number of Causal SNPs
    --effect  | -e    Effect size distribution.
                      0 for exponential
                      1 for Chi-Square
                      2 for Normal Distribution
    --fix     | -f    Fixed effect
    --extract | -E    List of SNPs to extract.
    --keep    | -k    List of samples to keep
    --thread  | -t    Number of thread used
    --herit   | -H    List of heritability to construct
    --seed    | -s    Seed for random generator
    --help    | -h    Display this help message
    
