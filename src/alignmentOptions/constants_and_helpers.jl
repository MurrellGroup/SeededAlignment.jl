# TODO replace with generation via bitmask
# NOTE important that they stay AminoAcid type to work with BioSequences.compatbits
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

## helper functions ##
# Convert NucleicAcid to integer: A -> 1, C -> 2, G -> 3, T -> 4
@inline toInt(x::NucleicAcid) = Int(trailing_zeros(BioSequences.compatbits(x)))+1 # decodes ambigious nucleotides badly
@inline toInt(aa::AminoAcid)::Int = Int(trailing_zeros(BioSequences.compatbits(aa)))+1
# translation function
@inline function fast_translate(dna_seq1::DNA,dna_seq2::DNA,dna_seq3::DNA)
    i1 = toInt(dna_seq1)
    i2 = toInt(dna_seq2)
    i3 = toInt(dna_seq3)
    hash_index = (i1 - 1)*16 + (i2 - 1)*4 + i3
    return SeededAlignment.CODON_TABLE[hash_index]
end

# test matrix

# Parameters
α = -0.4   # transition rate
β = -0.2   # transversion rate
# Explicit 4x4 K2P rate matrix Q
# Order: A C G T 
const NUC_SUB_MATRIX = [
    0   -β         -α       -β;
         -β   0    -β       -α;
         -α     -β     0    -β;
         -β     -α         -β   0
]