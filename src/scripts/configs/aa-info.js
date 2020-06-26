const STOP_CODONS = {
  UAA: {key: 'UAA', name: 'Stop codon UAA'},
  UAG: {key: 'UAG', name: 'Stop codon UAG'},
  UGA: {key: 'UGA', name: 'Stop codon UGA'},
};

export const AMINO_ACIDS = {
  X: {
    key: 'X', 
    name: 'Stop codon',
    nameAbbrev: 'Ter',
    codons: ['TAA','TAG','TGA'],
    color: '#333333', // dark grey
  },
  W: {
    key: 'W', 
    name: 'Tryptophan',
    nameAbbrev: 'Trp',
    essential: true,
    codons: ['TGG'],
    property: 'aromatic',
    color: '#abc745', // green
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/trp.png'
  },
  C: {
    key: 'C', 
    name: 'Cystine',
    nameAbbrev: 'Cys',
    essential: false,
    codons: ['TGT','TGC'],
    property: 'sulfur-containing',
    color: '#e9ae02', // orange
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/cys.png'
  },
  G: {
    key: 'G', 
    name: 'Glycine',
    nameAbbrev: 'Gly',
    essential: false,
    codons: ['GGT','GGC','GGA','GGG'],
    property: 'aliphatic',
    color: '#d81c24', // red
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/gly.png'
  },
  R: {
    key: 'R', 
    name: 'Arginine',
    nameAbbrev: 'Arg',
    essential: false,
    codons: ['CGT','CGC','CGA','CGG','AGA','AGG'],
    property: 'basic',
    color: '#739ed4', // light blue
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/arg.png'
  },
  S: {
    key: 'S', 
    name: 'Serine',
    nameAbbrev: 'Ser',
    essential: false,
    codons: ['TCT','TCC','TCA','TCG','AGT','AGC'],
    property: 'hydroxylic',
    color: '#f39bbf', // pink
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/ser.png'
  },
  T: {
    key: 'T', 
    name: 'Threonine',
    nameAbbrev: 'Thr',
    essential: true,
    codons: ['ACT','ACC','ACA','ACG'],
    property: 'hydroxylic',
    color: '#f39bbf', // pink
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/thr.png'
  },
  A: {
    key: 'A', 
    name: 'Alanine',
    nameAbbrev: 'Ala',
    essential: false,
    codons: ['GCT','GCC','GCA','GCG'],
    property: 'aliphatic',
    color: '#d81c24', // red
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/ala.png'
  },
  P: {
    key: 'P', 
    name: 'Proline',
    nameAbbrev: 'Pro',
    essential: false,
    codons: ['CCT','CCC','CCA','CCG'],
    property: 'aliphatic',
    color: '#d81c24', // red
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/pro.png'
  },
  F: {
    key: 'F', 
    name: 'Phenylalanine',
    nameAbbrev: 'Phe',
    essential: true,
    codons: ['TTT','TTC'],
    property: 'aromatic',
    color: '#abc745', // green
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/phe.png'
  },
  L: {
    key: 'L', 
    name: 'Leucine',
    nameAbbrev: 'Leu',
    essential: true,
    codons: ['CTT','CTC','CTA','CTG','TTA','TTG'],
    property: 'aliphatic',
    color: '#d81c24', // red
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/leu.png'
  },
  V: {
    key: 'V', 
    name: 'Valine',
    nameAbbrev: 'Val',
    essential: true,
    codons: ['GTT','GTC','GTA','GTG'],
    property: 'aliphatic',
    color: '#d81c24', // red
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/val.png'
  },
  I: {
    key: 'I', 
    name: 'Isoleucine',
    nameAbbrev: 'Ile',
    essential: true,
    codons: ['ATT','ATC','ATA'],
    property: 'aliphatic',
    color: '#d81c24', // red
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/ile.png'
  },
  M: {
    key: 'M', 
    name: 'Methionine',
    nameAbbrev: 'Met',
    essential: true,
    codons: ['ATG'],
    property: 'sulfur-containing',
    color: '#e9ae02', // orange
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/met.png'
  },
  Q: {
    key: 'Q', 
    name: 'Glutamine',
    nameAbbrev: 'Gln',
    essential: false,
    codons: ['CAA','CAG'],
    property: 'amidic',
    color: '#2c3385', // dark blue
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/gln.png'
  },
  H: {
    key: 'H', 
    name: 'Histidine',
    nameAbbrev: 'His',
    essential: true,
    codons: ['CAT','CAC'],
    property: 'basic',
    color: '#739ed4', // light blue
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/his.png'
  },
  D: {
    key: 'D', 
    name: 'Aspartic acid',
    nameAbbrev: 'Asp',
    essential: false,
    codons: ['GAT','GAC'],
    property: 'acidic',
    color: '#e96b17', // dark orange
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/asp.png'
  },
  E: {
    key: 'E', 
    name: 'Glutamic acid',
    nameAbbrev: 'Glu',
    essential: false,
    codons: ['GAA','GAG'],
    property: 'acidic',
    color: '#e96b17', // dark orange
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/glu.png'
  },
  K: {
    key: 'K', 
    name: 'Lysine',
    nameAbbrev: 'Lys',
    essential: true,
    codons: ['AAA','AAG'],
    property: 'basic',
    color: '#739ed4', // light blue
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/lys.png'
  },
  N: {
    key: 'N', 
    name: 'Asparagine',
    nameAbbrev: 'Asn',
    essential: false,
    codons: ['AAT','AAC'],
    property: 'amidic',
    color: '#2c3385', // dark blue
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/asn.png'
  },
  Y: {
    key: 'Y', 
    name: 'Tyrosine',
    nameAbbrev: 'Tyr',
    essential: false,
    codons: ['TAT','TAC'],
    property: 'aromatic',
    color: '#abc745', // green
    image: 'https://aveit.s3.amazonaws.com/higlass/static/amino-acids/tyr.png'
  },
};

export const CODONS = {
  TTT: AMINO_ACIDS.F,
  TTC: AMINO_ACIDS.F,
  TTA: AMINO_ACIDS.L,
  TTG: AMINO_ACIDS.L,
  TCT: AMINO_ACIDS.S,
  TCC: AMINO_ACIDS.S,
  TCA: AMINO_ACIDS.S,
  TCG: AMINO_ACIDS.S,
  TAT: AMINO_ACIDS.Y,
  TAC: AMINO_ACIDS.Y,
  TAA: AMINO_ACIDS.X,
  TAG: AMINO_ACIDS.X,
  TGT: AMINO_ACIDS.C,
  TGC: AMINO_ACIDS.C,
  TGA: AMINO_ACIDS.X,
  TGG: AMINO_ACIDS.W,
  CTT: AMINO_ACIDS.L,
  CTC: AMINO_ACIDS.L,
  CTA: AMINO_ACIDS.L,
  CTG: AMINO_ACIDS.L,
  CCT: AMINO_ACIDS.P,
  CCC: AMINO_ACIDS.P,
  CCA: AMINO_ACIDS.P,
  CCG: AMINO_ACIDS.P,
  CAT: AMINO_ACIDS.H,
  CAC: AMINO_ACIDS.H,
  CAA: AMINO_ACIDS.Q,
  CAG: AMINO_ACIDS.Q,
  CGT: AMINO_ACIDS.R,
  CGC: AMINO_ACIDS.R,
  CGA: AMINO_ACIDS.R,
  CGG: AMINO_ACIDS.R,
  ATT: AMINO_ACIDS.I,
  ATC: AMINO_ACIDS.I,
  ATA: AMINO_ACIDS.I,
  ATG: AMINO_ACIDS.M,
  ACT: AMINO_ACIDS.T,
  ACC: AMINO_ACIDS.T,
  ACA: AMINO_ACIDS.T,
  ACG: AMINO_ACIDS.T,
  AAT: AMINO_ACIDS.N,
  AAC: AMINO_ACIDS.N,
  AAA: AMINO_ACIDS.K,
  AAG: AMINO_ACIDS.K,
  AGT: AMINO_ACIDS.S,
  AGC: AMINO_ACIDS.S,
  AGA: AMINO_ACIDS.R,
  AGG: AMINO_ACIDS.R,
  GTT: AMINO_ACIDS.V,
  GTC: AMINO_ACIDS.V,
  GTA: AMINO_ACIDS.V,
  GTG: AMINO_ACIDS.V,
  GCT: AMINO_ACIDS.A,
  GCC: AMINO_ACIDS.A,
  GCA: AMINO_ACIDS.A,
  GCG: AMINO_ACIDS.A,
  GAT: AMINO_ACIDS.D,
  GAC: AMINO_ACIDS.D,
  GAA: AMINO_ACIDS.E,
  GAG: AMINO_ACIDS.E,
  GGT: AMINO_ACIDS.G,
  GGC: AMINO_ACIDS.G,
  GGA: AMINO_ACIDS.G,
  GGG: AMINO_ACIDS.G,
};

export default AMINO_ACIDS;
