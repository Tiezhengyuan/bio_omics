'''
https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/protein-biology/protein-structural-analysis/amino-acid-reference-chart#hydrophobicity

ApH 2 values: Normalized from Sereda et al., J. Chrom. 676: 139-153 (1994).
BpH 7 values: Monera et al., J. Protein Sci. 1: 319-329 (1995).

https://pmc.ncbi.nlm.nih.gov/articles/PMC5037676/
'''
# hydrophobicity index at pH7
AA = {
    "A":{'abb': 'ala', 'hydrophobicity': 41,},
    "R":{'abb': 'arg', 'hydrophobicity': -14,},
    "N":{'abb': 'asn', 'hydrophobicity': -28,},
    "D":{'abb': 'asp', 'hydrophobicity': -55,},
    "C":{'abb': 'cys', 'hydrophobicity': 49,},
    "E":{'abb': 'glu', 'hydrophobicity': -31,},
    "Q":{'abb': 'gln', 'hydrophobicity': -10,},
    "G":{'abb': 'gly', 'hydrophobicity': 0,},
    "H":{'abb': 'his', 'hydrophobicity': 8,},
    "I":{'abb': 'ile', 'hydrophobicity': 99,},
    "L":{'abb': 'leu', 'hydrophobicity': 97,},
    "K":{'abb': 'lys', 'hydrophobicity': -23,},
    "M":{'abb': 'met', 'hydrophobicity': 74,},
    "F":{'abb': 'phe', 'hydrophobicity': 100,},
    "P":{'abb': 'pro', 'hydrophobicity': -46,}, #PH2
    "S":{'abb': 'ser', 'hydrophobicity': -5,},
    "T":{'abb': 'thr', 'hydrophobicity': 13,},
    "W":{'abb': 'trp', 'hydrophobicity': 97,},
    "Y":{'abb': 'tyr', 'hydrophobicity': 63,},
    "V":{'abb': 'val', 'hydrophobicity': 76,},
}

PROPERTY = {
    'A': {
        'name'                  : 'Alanine',
        'abbreviation'          : 'Ala',
        'hydrophobicity_ph7'    : 41,
        'hydrophobicity'        : 0.61,
        'polarity'              : 8.1,
        'polarizability'        : 0.046,
        'van_der_Waals_volume'  : 1,
    },
    'R': {
        'name'                  : 'Arginine',
        'abbreviation'          : 'Arg',
        'hydrophobicity_ph7'    : -14,
        'hydrophobicity'        : 0.60,
        'polarity'              : 10.5,
        'polarizability'        : 0.291,
        'van_der_Waals_volume'  : 6.13,
    },
    'N': {
        'name'                  : 'Asparagine',
        'abbreviation'          : 'Asn',
        'hydrophobicity_ph7'    : -28,
        'hydrophobicity'        : 0.06,
        'polarity'              : 11.6,
        'polarizability'        : 0.134,
        'van_der_Waals_volume'  : 2.95,
    },
    'D': {
        'name'                  : 'Aspartic Acid',
        'abbreviation'          : 'Asp',
        'hydrophobicity_ph7'    : -55,
        'hydrophobicity'        : 0.46,
        'polarity'              : 13.0,
        'polarizability'        : 0.105,
        'van_der_Waals_volume'  : 2.78,
    },
    'C': {'name'                : 'Cysteine',
        'abbreviation'          : 'Cys',
        'hydrophobicity_ph7'    : 49,
        'hydrophobicity'        : 1.07,
        'polarity'              : 5.5,
        'polarizability'        : 0.128,
        'van_der_Waals_volume'  : 2.43,
    },
    'Q': {
        'name'                  : 'Glutamine',
        'abbreviation'          : 'Gln',
        'hydrophobicity_ph7'    : -10,
        'hydrophobicity'        : 0,
        'polarity'              : 10.5,
        'polarizability'        : 0.180,
        'van_der_Waals_volume'  : 3.95,
    },
    'E': {
        'name'                  : 'Glutamic Acid',
        'abbreviation'          : 'Glu',
        'hydrophobicity_ph7'    : -31,
        'hydrophobicity'        : 0.47,
        'polarity'              : 12.3,
        'polarizability'        : 0.151,
        'van_der_Waals_volume'  : 3.78,
    },
    'G': {
        'name'                  : 'Glycine',
        'abbreviation'          : 'Gly',
        'hydrophobicity_ph7'    : 0,
        'hydrophobicity'        : 0.07,
        'polarity'              : 9.0,
        'polarizability'        : 0,
        'van_der_Waals_volume'  : 0,
    },
    'H': {
        'name'                  : 'Histidine',
        'abbreviation'          : 'His',
        'hydrophobicity_ph7'    : 8,
        'hydrophobicity'        : 0.61,
        'polarity'              : 10.4,
        'polarizability'        : 0.23,
        'van_der_Waals_volume'  : 4.66
    },
    'I': {
        'name'                  : 'Isoleucine',
        'abbreviation'          : 'Ile',
        'hydrophobicity_ph7'    : 99,
        'hydrophobicity'        : 2.22,
        'polarity'              : 5.2,
        'polarizability'        : 0.186,
        'van_der_Waals_volume'  : 4,
    },
    'L': {
        'name'                  : 'Leucine',
        'abbreviation'          : 'Leu',
        'hydrophobicity_ph7'    : 97,
        'hydrophobicity'        : 1.53,
        'polarity'              : 4.9,
        'polarizability'        : 0.186,
        'van_der_Waals_volume'  : 4,
    },
    'K': {
        'name'                  : 'Lysine',
        'abbreviation'          : 'Lys',
        'hydrophobicity_ph7'    : -23,
        'hydrophobicity'        : 1.15,
        'polarity'              : 11.3,
        'polarizability'        : 0.219,
        'van_der_Waals_volume'  : 4.77,
    },
    'M': {
        'name'                  : 'Methionine',
        'abbreviation'          : 'Met',
        'hydrophobicity_ph7'    : 74,
        'hydrophobicity'        : 1.18,
        'polarity'              : 5.7,
        'polarizability'        : 0.221,
        'van_der_Waals_volume'  : 4.43
    },
    'F': {
        'name'                  : 'Phenylalanine',
        'abbreviation'          : 'Phe',
        'hydrophobicity_ph7'    : 100,
        'hydrophobicity'        : 2.02,
        'polarity'              : 5.2,
        'polarizability'        : 0.29,
        'van_der_Waals_volume'  : 5.89,
    },
    'P': {
        'name'                  : 'Proline',
        'abbreviation'          : 'Pro',
        'hydrophobicity_ph7'    : -46,
        'hydrophobicity'        : 1.95,
        'polarity'              : 8,
        'polarizability'        : 0.131,
        'van_der_Waals_volume'  : 2.72,
    },
    'S': {
        'name'                  : 'Serine',
        'abbreviation'          : 'Ser',
        'hydrophobicity_ph7'    : -5,
        'hydrophobicity'        : 0.05,
        'polarity'              : 9.2,
        'polarizability'        : 0.062,
        'van_der_Waals_volume'  : 1.6,
    },
    'T': {
        'name'                  : 'Threonine',
        'abbreviation'          : 'Thr',
        'hydrophobicity_ph7'    : 13,
        'hydrophobicity'        : 0.05,
        'polarity'              : 8.6,
        'polarizability'        : 0.108,
        'van_der_Waals_volume'  : 2.6,
    },
    'W': {
        'name'                  : 'Tryptophan',
        'abbreviation'          : 'Trp',
        'hydrophobicity_ph7'    : 97,
        'hydrophobicity'        : 2.65,
        'polarity'              : 5.4,
        'polarizability'        : 0.409,
        'van_der_Waals_volume'  : 8.08,
    },
    'Y': {
        'name'                  : 'Tyrosine',
        'abbreviation'          : 'Tyr',
        'hydrophobicity_ph7'    : 63,
        'hydrophobicity'        : 1.88,
        'polarity'              : 6.2,
        'polarizability'        : 0.298,
        'van_der_Waals_volume'  : 6.47,
    },
    'V': {
        'name'                  : 'Valine',
        'abbreviation'          : 'Val',
        'hydrophobicity_ph7'    : 76,
        'hydrophobicity'        : 1.32,
        'polarity'              : 5.9,
        'polarizability'        : 0.140,
        'van_der_Waals_volume'  : 3,
    },
}

SMILES = {
    "A": 'CC(C(=O)O)N',
    "R": 'C(CC(=N)N)CC(C(=O)O)N',
    "N": 'C(C(C(=O)O)N)C(=O)N',
    "D": 'C(C(C(=O)O)N)C(=O)O',
    "C": 'C(C(C(=O)O)N)S',
    "E": 'C(CC(=O)O)C(C(=O)O)N',
    "Q": 'C(CC(=O)N)C(C(=O)O)N',
    'G': 'C(C(=O)O)N',
    "H": 'C1=C(NC=N1)CC(C(=O)O)N',
    "I": 'CCC(C)C(C(=O)O)N',
    "L": 'CC(C)CC(C(=O)O)N',
    "K": 'C(CCN)CC(C(=O)O)N',    
    "M": 'CSCC(C(=O)O)N',
    "F": 'C1=CC=C(C=C1)CC(C(=O)O)N',
    "P": 'C1CC(NC1)C(=O)O',
    "S": 'C(C(C(=O)O)N)O',
    "T": 'CC(C(C(=O)O)N)O',
    "W": 'C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N',
    'Y': 'C1=CC(=CC=C1CC(C(=O)O)N)O',
    "V": 'CC(C)C(C(=O)O)N',
}
'''
Linear formula
Alanine	Ala	A	C3H7NO2	CH3-CH(NH2)-COOH
Arginine	Arg	R	C6H14N4O2	HN=C(NH2)-NH-(CH2)3-CH(NH2)-COOH
Asparagine	Asn	N	C4H8N2O3	H2N-CO-CH2-CH(NH2)-COOH
Aspartic acid	Asp	D	C4H7NO4	HOOC-CH2-CH(NH2)-COOH
Cysteine	Cys	C	C3H7NO2S	HS-CH2-CH(NH2)-COOH
Glutamine	Gln	Q	C5H10N2O3	H2N-CO-(CH2)2-CH(NH2)-COOH
Glutamic acid	Glu	E	C5H9NO4	HOOC-(CH2)2-CH(NH2)-COOH
Glycine	Gly	G	C2H5NO2	NH2-CH2-COOH
Histidine	His	H	C6H9N3O2	NH-CH=N-CH=C-CH2-CH(NH2)-COOH
Isoleucine	Ile	I	C6H13NO2	CH3-CH2-CH(CH3)-CH(NH2)-COOH
Leucine	Leu	L	C6H13NO2	(CH3)2-CH-CH2-CH(NH2)-COOH
Lysine	Lys	K	C6H14N2O2	H2N-(CH2)4-CH(NH2)-COOH
Methionine	Met	M	C5H11NO2S	CH3-S-(CH2)2-CH(NH2)-COOH
Phenylalanine	Phe	F	C9H11NO2	Ph-CH2-CH(NH2)-COOH
Proline	Pro	P	C5H9NO2	NH-(CH2)3-CH-COOH
Serine	Ser	S	C3H7NO3	HO-CH2-CH(NH2)-COOH
Threonine	Thr	T	C4H9NO3	CH3-CH(OH)-CH(NH2)-COOH
Tryptophan	Trp	W	C11H12N2O2	Ph-NH-CH=C-CH2-CH(NH2)-COOH
Tyrosine	Tyr	Y	C9H11NO3	HO-Ph-CH2-CH(NH2)-COOH
Valine	Val	V	C5H11NO2	(CH3)2-CH-CH(NH2)-COOH
'''
SMILES_1 = {
    'A': 'C',
    'R': 'MBBB',
    'N': 'NJB',
    'D': 'EB',
    'C': 'KB',
    'Q': 'NJBB',
    'E': 'EBB',
    'G': '',
    'H': 'HB',
    'I': 'CBI',
    'L': 'DAB',
    'K': 'NBBBB',
    'M': 'CSBB',
    'F': 'FB',
    'P': 'P',
    'S': 'S',
    'T': 'CT',
    'W': 'FWB',
    'Y': 'OFB',
    'V': 'DA',
}
