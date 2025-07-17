# useful constants
const AA_STANDARD_ORDER = AminoAcid.([
    'A','R','N','D','C','Q','E','G','H','I',
    'L','K','M','F','P','S','T','W','Y','V',
    '*'
])

const CODON_TABLE::Vector{AminoAcid} = AminoAcid[
    # AAA AAC AAG AAT
    AminoAcid('K'), AminoAcid('N'), AminoAcid('K'), AminoAcid('N'),
    # ACA ACC ACG ACT
    AminoAcid('T'), AminoAcid('T'), AminoAcid('T'), AminoAcid('T'),
    # AGA AGC AGG AGT
    AminoAcid('R'), AminoAcid('S'), AminoAcid('R'), AminoAcid('S'),
    # ATA ATC ATG ATT
    AminoAcid('I'), AminoAcid('I'), AminoAcid('M'), AminoAcid('I'),

    # CAA CAC CAG CAT
    AminoAcid('Q'), AminoAcid('H'), AminoAcid('Q'), AminoAcid('H'),
    # CCA CCC CCG CCT
    AminoAcid('P'), AminoAcid('P'), AminoAcid('P'), AminoAcid('P'),
    # CGA CGC CGG CGT
    AminoAcid('R'), AminoAcid('R'), AminoAcid('R'), AminoAcid('R'),
    # CTA CTC CTG CTT
    AminoAcid('L'), AminoAcid('L'), AminoAcid('L'), AminoAcid('L'),

    # GAA GAC GAG GAT
    AminoAcid('E'), AminoAcid('D'), AminoAcid('E'), AminoAcid('D'),
    # GCA GCC GCG GCT
    AminoAcid('A'), AminoAcid('A'), AminoAcid('A'), AminoAcid('A'),
    # GGA GGC GGG GGT
    AminoAcid('G'), AminoAcid('G'), AminoAcid('G'), AminoAcid('G'),
    # GTA GTC GTG GTT
    AminoAcid('V'), AminoAcid('V'), AminoAcid('V'), AminoAcid('V'),

    # TAA TAC TAG TAT
    AminoAcid('*'), AminoAcid('Y'), AminoAcid('*'), AminoAcid('Y'),
    # TCA TCC TCG TCT
    AminoAcid('S'), AminoAcid('S'), AminoAcid('S'), AminoAcid('S'),
    # TGA TGC TGG TGT
    AminoAcid('*'), AminoAcid('C'), AminoAcid('W'), AminoAcid('C'),
    # TTA TTC TTG TTT
    AminoAcid('L'), AminoAcid('F'), AminoAcid('L'), AminoAcid('F')
]

# score_matrix for the 20 amino acids - so no stop codon 
# TODO add stop_codon scoring
const BLOSUM62::Matrix{Float64} = [
    4.0 -1.0 -2.0 -2.0  0.0 -1.0 -1.0  0.0 -2.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0  1.0  0.0 -3.0 -2.0  0.0;
   -1.0  5.0  0.0 -2.0 -3.0  1.0  0.0 -2.0  0.0 -3.0 -2.0  2.0 -1.0 -3.0 -2.0 -1.0 -1.0 -3.0 -2.0 -3.0;
   -2.0  0.0  6.0  1.0 -3.0  0.0  0.0  0.0  1.0 -3.0 -3.0  0.0 -2.0 -3.0 -2.0  1.0  0.0 -4.0 -2.0 -3.0;
   -2.0 -2.0  1.0  6.0 -3.0  0.0  2.0 -1.0 -1.0 -3.0 -4.0 -1.0 -3.0 -3.0 -1.0  0.0 -1.0 -4.0 -3.0 -3.0;
    0.0 -3.0 -3.0 -3.0  9.0 -3.0 -4.0 -3.0 -3.0 -1.0 -1.0 -3.0 -1.0 -2.0 -3.0 -1.0 -1.0 -2.0 -2.0 -1.0;
   -1.0  1.0  0.0  0.0 -3.0  5.0  2.0 -2.0  0.0 -3.0 -2.0  1.0  0.0 -3.0 -1.0  0.0 -1.0 -2.0 -1.0 -2.0;
   -1.0  0.0  0.0  2.0 -4.0  2.0  5.0 -2.0  0.0 -3.0 -3.0  1.0 -2.0 -3.0 -1.0  0.0 -1.0 -3.0 -2.0 -2.0;
    0.0 -2.0  0.0 -1.0 -3.0 -2.0 -2.0  6.0 -2.0 -4.0 -4.0 -2.0 -3.0 -3.0 -2.0  0.0 -2.0 -2.0 -3.0 -3.0;
   -2.0  0.0  1.0 -1.0 -3.0  0.0  0.0 -2.0  8.0 -3.0 -3.0 -1.0 -2.0 -1.0 -2.0 -1.0 -2.0 -2.0  2.0 -3.0;
   -1.0 -3.0 -3.0 -3.0 -1.0 -3.0 -3.0 -4.0 -3.0  4.0  2.0 -3.0  1.0  0.0 -3.0 -2.0 -1.0 -3.0 -1.0  3.0;
   -1.0 -2.0 -3.0 -4.0 -1.0 -2.0 -3.0 -4.0 -3.0  2.0  4.0 -2.0  2.0  0.0 -3.0 -2.0 -1.0 -2.0 -1.0  1.0;
   -1.0  2.0  0.0 -1.0 -3.0  1.0  1.0 -2.0 -1.0 -3.0 -2.0  5.0 -1.0 -3.0 -1.0  0.0 -1.0 -3.0 -2.0 -2.0;
   -1.0 -1.0 -2.0 -3.0 -1.0  0.0 -2.0 -3.0 -2.0  1.0  2.0 -1.0  5.0  0.0 -2.0 -1.0 -1.0 -1.0 -1.0  1.0;
   -2.0 -3.0 -3.0 -3.0 -2.0 -3.0 -3.0 -3.0 -1.0  0.0  0.0 -3.0  0.0  6.0 -4.0 -2.0 -2.0  1.0  3.0 -1.0;
   -1.0 -2.0 -2.0 -1.0 -3.0 -1.0 -1.0 -2.0 -2.0 -3.0 -3.0 -1.0 -2.0 -4.0  7.0 -1.0 -1.0 -4.0 -3.0 -2.0;
    1.0 -1.0  1.0  0.0 -1.0  0.0  0.0  0.0 -1.0 -2.0 -2.0  0.0 -1.0 -2.0 -1.0  4.0  1.0 -3.0 -2.0 -2.0;
    0.0 -1.0  0.0 -1.0 -1.0 -1.0 -1.0 -2.0 -2.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0  1.0  5.0 -2.0 -2.0  0.0;
   -3.0 -3.0 -4.0 -4.0 -2.0 -2.0 -3.0 -2.0 -2.0 -3.0 -2.0 -3.0 -1.0  1.0 -4.0 -3.0 -2.0 11.0  2.0 -3.0;
   -2.0 -2.0 -2.0 -3.0 -2.0 -1.0 -2.0 -3.0  2.0 -1.0 -1.0 -2.0 -1.0  3.0 -3.0 -2.0 -2.0  2.0  7.0 -1.0;
    0.0 -3.0 -3.0 -3.0 -1.0 -2.0 -2.0 -3.0 -3.0  3.0  1.0 -2.0  1.0 -1.0 -2.0 -2.0  0.0 -3.0 -1.0  4.0
]

# TODO add ambigious nucleotide scoring
# Order: A, C, G, T
const NUC_MATRIX = [
     1.0  -2.0   0.0  -2.0;  # A
    -2.0   1.0  -2.0   0.0;  # C
     0.0  -2.0   1.0  -2.0;  # G
    -2.0   0.0  -2.0   1.0   # T
]

## helper functions ##

# TODO unit test toInt functions
# Convert NucleicAcid to integer: A -> 1, C -> 2, G -> 3, T -> 4
@inline toInt(x::NucleicAcid) = Int(trailing_zeros(BioSequences.compatbits(x)))+1 # decodes ambigious nucleotides badly
@inline toInt(aa::AminoAcid)::Int = Int(trailing_zeros(BioSequences.compatbits(aa)))+1
# TODO unit test the CODON_TABLE
@inline function fast_translate(dna_seq::NTuple{3, DNA})
    i1 = toInt(dna_seq[1])
    i2 = toInt(dna_seq[2])
    i3 = toInt(dna_seq[3])
    hash_index = (i1 - 1)*16 + (i2 - 1)*4 + i3
    return SeededAlignment.CODON_TABLE[hash_index]
end