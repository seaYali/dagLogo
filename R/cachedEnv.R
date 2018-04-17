## create a cached environment to store those package-wide accessible constants.
cachedEnv <- new.env(parent = emptyenv())

## color and character symbol encodings and grouping for amino acids
no <- list(
    color = c('A' = '#e6194b', 'C' = '#3cb44b', 'D' = '#0082c8',
              'E' = '#f58231', 'F' = '#911eb4', 'G' = '#46f0f0',
              'H' = '#f032e6', 'I' = '#d2f53c', 'K' = '#fabebe',
              'L' = '#008080', 'M' = '#e6beff', 'N' = '#aa6e28',
              'P' = '#00FF00', 'Q' = '#800000', 'R' = '#aaffc3',
              'S' = '#808000', 'T' = '#ffd8b1', 'V' = '#000080',
              'W' = '#808080', 'Y' = '#000000'),
    symbol = c("Ala" = "A", "Arg" = "R", "Asn" = "N", "Asp" = "D",
               "Cys" = "C", "Glu" = "E", "Gln" = "Q", "Gly" = "G",
               "His" = "H", "Ile" = "I", "Leu" = "L", "Lys" = "K",
               "Met" = "M", "Phe" = "F", "Pro" = "P", "Ser" = "S",
               "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V"),
    group = NULL)

classic <- list(
    color = c(
        nonpolar_aliphatic = "#000000",
        polar_uncharged = "#00811B",
        aromatic = "#2000C7",
        positively_charged = "#800080",
        negatively_charged = "#FFB32C"),
    symbol = c(
        nonpolar_aliphatic = "M",
        polar_uncharged = "U",
        aromatic = "A",
        positively_charged = "P",
        negatively_charged = "N"),
    ## Aliphatic: pertaining to nonaromatic hydrocarbon compounds in which the 
    ## constituent carbon atoms can be straight-chain, branched chain, or 
    ## cyclic, as in alicyclic compounds; saturated, as in the paraffins; or 
    ## unsaturated, as in the olefins and alkynes. (Ref. Betts and Russell. 
    ## Amino acid properties and consequences of substitutions.  
    ## Bioinformtics for Geneticists. 2003)
    group = list(
        nonpolar_aliphatic = c("A", "G", "I", "L", "M","V","P"),
        polar_uncharged = c("C","S", "T", "N", "Q"),
        aromatic = c("F", "W", "Y"),
        positively_charged = c("H", "K", "R"),
        negatively_charged = c("D", "E"))
)

chemistry <- list(
    color = c(
        hydrophobic = "#000000",
        polar = "#00811B",
        basic = "#2000C7",
        neutral  = "#800080",
        acidic = "#D00001"),
    symbol = c(
        hydrophobic = "H",
        polar = "P",
        basic = "B",
        neutral  = "U",
        acidic = "A"),
    group = list(
        hydrophobic = c("G", "A", "I", "L", "V", "M", "P", "F","W"),
        polar = c("C", "S", "T", "Y"),
        basic = c("H", "K", "R"),
        neutral = c("N", "Q"),
        acidic = c("D", "E"))
)

hydrophobicity <- list(
    color = c(
        hydrophilic = '#000000',
        neutral = '#00811B',
        hydrophobic  = '#2000C7'),
    symbol = c(
        hydrophilic = 'I',
        neutral = 'U',
        hydrophobic = 'O'),
    group = list(
        hydrophilic = c("D", "E", "K", "R", "N", "Q"),
        neutral = c("A", "G", "H", "P", "S", "T"),
        hydrophobic = c("C", "F", "I", "L", "M", "V", "W", "Y"))
)

charge <- list(
    color = c(
        positive = "#FFB32C",
        neutral = "#2000C7",
        negative = "#CCCCCC"),
    symbol = c(
        positive = "P",
        neutral = "U",
        negative = "N"),
    group = list(
        positive = c("H", "K", "R"),
        neutral = c("A", "C", "F", "G", "I", "L", "M", "N", "P", 
                    "Q", "S", "T", "V","W","Y"),
        negative = c("D", "E"))
) 

## 8-level schemes for amino acid alphabet reduction derived from correlations 
## based on the BLOSUM50 similarity matrix. More schemes are available from
## Stephenson and Freeland. Unwearthing the root of amino acid similarity. 
## J Mol. Evol. 2013 (77):159-169. Smith and Smith. Protein Eng. 1992(5):35-41.
## And many other literature.

BLOSM50_L1 <- list(
    color = c(
        LVIMCAGSTPFYW = "#008080",
        EDNQKRH = "#f58231"),
    symbol = c(
        LVIMCAGSTPFYW = "L",
        EDNQKRH = "E"),
    group = list(
        LVIMCAGSTPFYW = c("L", "V", "I", "M", "C", "A", "G", "S", "T", "P", 
                          "F", "Y", "W"),
        EDNQKRH = c("E", "D", "N", "Q", "K", "R", "H"))
) 

BLOSM50_L2 <- list(
    color = c(
        LVIMCAGSTP = "#008080",
        FYW = "#911eb4",
        EDNQKRH = "#f58231"),
    symbol = c(
        LVIMCAGSTP = "L",
        FYW = "F",
        EDNQKRH = "E"),
    group = list(
        LVIMCAGSTP = c("L", "V", "I", "M", "C", "A", "G", "S", "T", "P"), 
        FYW = c("F", "Y", "W"),
        EDNQKRH = c("E", "D", "N", "Q", "K", "R", "H"))
) 

BLOSM50_L3 <- list(
    color = c(
        LVIMC = "#008080",
        AGSTP = "#e6194b",
        FYW = '#911eb4',
        EDNQKRH = "#f58231"),
    symbol = c(
        LVIMC = "L",
        AGSTP = "A",
        FYW = "F",
        EDNQKRH = "E"),
    group = list(
        LVIMC = c("L", "V", "I", "M", "C"), 
        AGSTP = c("A", "G", "S", "T", "P"),
        FYW = c("F", "Y", "W"),
        EDNQKRH = c("E", "D", "N", "Q", "K", "R", "H"))
) 

BLOSM50_L4 <- list(
    color = c(
        LVIMC = "#008080",
        AGSTP = "#e6194b",
        FYW = '#911eb4',
        EDNQ = "#f58231",
        KRH = "#fabebe"),
    symbol = c(
        LVIMC = "L",
        AGSTP = "A",
        FYW = "F",
        EDNQ = "E",
        KRH = "K"),
    group = list(
        LVIMC = c("L", "V", "I", "M", "C"), 
        AGSTP = c("A", "G", "S", "T", "P"),
        FYW = c("F", "Y", "W"),
        EDNQ = c("E", "D", "N", "Q"),
        KRH = c("K", "R", "H"))
) 

BLOSM50_L5 <- list(
    color = c(
        LVIMC = "#008080",
        AGST = "#e6194b",
        P = "#00FF00",
        FYW = "#911eb4",
        EDNQ = "#f58231",
        KRH = "#fabebe"),
    symbol = c(
        LVIMC = "L",
        AGST = "A",
        P = "P",
        FYW = "F",
        EDNQ = "E",
        KRH = "K"),
    group = list(
        LVIMC = c("L", "V", "I", "M", "C"), 
        AGST = c("A", "G", "S", "T"),
        P = "P",
        FYW = c("F", "Y", "W"),
        EDNQ = c("E", "D", "N", "Q"),
        KRH = c("K", "R", "H"))
) 


BLOSM50_L6 <- list(
    color = c(
        LVIMC = "#008080",
        AG = "#e6194b",
        ST = "#808000",
        P = "#00FF00",
        FYW = "#911eb4",
        EDNQ = "#f58231",
        KR = "#fabebe",
        H = "#f032e6"),
    symbol = c(
        LVIMC = "L",
        AG = "A",
        ST = "S",
        P = "P",
        FYW = "F",
        EDNQ = "E",
        KR = "K",
        H = "H"),
    group = list(
        LVIMC = c("L", "V", "I", "M", "C"), 
        AG = c("A", "G"),
        ST = c("S", "T"),
        P = "P",
        FYW = c("F", "Y", "W"),
        EDNQ = c("E", "D", "N", "Q"),
        KR = c("K", "R"),
        H = "H")
) 

BLOSM50_L7 <- list(
    color = c(
        LVIM = "#008080",
        C = "#3cb44b",
        A = "#e6194b",
        G = "#46f0f0",
        ST = "#808000",
        P = "#00FF00",
        FYW = "#911eb4",
        EDNQ = "#f58231",
        KR = "#fabebe",
        H = "#f032e6"),
    symbol = c(
        LVIM = "L",
        C = "C",
        A = "A",
        G = "G",
        ST = "S",
        P = "P",
        FYW = "F",
        EDNQ = "E",
        KR = "K",
        H = "H"),
    group = list(
        LVIM = c("L", "V", "I", "M"), 
        C = "C",
        A = "A",
        G = "G",
        ST = c("S", "T"),
        P = "P",
        FYW = c("F", "Y", "W"),
        EDNQ = c("E", "D", "N", "Q"),
        KR = c("K", "R"),
        H = "H")
) 

BLOSM50_L8 <- list(
    color = c(
        LVIM = "#008080",
        C = "#3cb44b",
        A = "#e6194b",
        G = "#46f0f0",
        S = "#808000",
        T = "#ffd8b1",
        P = "#00FF00",
        FY = "#911eb4",
        W = "#808080",
        E = "#f58231",
        D = "#0082c8",
        N = "#aa6e28",
        Q = "#800000",
        KR = "#fabebe",
        H = "#f032e6"),
    symbol = c(
        LVIM = "L",
        C = "C",
        A = "A",
        G = "G",
        S = "S",
        T = "T",
        P = "P",
        FY = "F",
        W = "W",
        E = "E",
        D = "D",
        N = "N",
        Q = "Q",
        KR = "K",
        H = "H"),
    group = list(
        LVIM = c("L", "V", "I", "M"), 
        C = "C",
        A = "A",
        G = "G",
        S = "S",
        T = "T",
        P = "P",
        FY = c("F", "Y"),
        W = "W",
        E = "E",
        D = "D",
        N = "N",
        Q = "Q",
        KR = c("K", "R"),
        H = "H")
) 


## Spatial frequency
contact_potential_Maiorov <- list(
    color = c(
        DENQ = "#0082c8",
        KR = "#fabebe",
        G = '#46f0f0',
        P = "#00FF00",
        AV = "#e6194b",
        STHWY = "#808000",
        CFMLI = "#3cb44b"),
    symbol = c(
        DENQ = "D",
        KR = "K",
        G = "G",
        P = "P",
        AV = "A",
        STHWY = "S",
        CFMLI = "C"),
    group = list(
        DENQ = c("D", "E", "N", "Q"), 
        KR = c("K", "R"),
        G = "G",
        P = "P",
        AV = c("A", "V"),
        STHWY = c("S", "T", "H", "W", "Y"),
        CFMLI = c("C", "F", "M", "L", "I"))
) 

## protein blocks
protein_blocks_Rogov <- list(
    color = c(
        DNSTA = "#0082c8",
        EKRQ = "#f58231",
        G = '#46f0f0',
        P = "#00FF00",
        H = "#f032e6",
        C = "#3cb44b",
        W = "#808080",
        M = "#e6beff",
        YFLIV = "#000000"),
    symbol = c(
        DNSTA = "D",
        EKRQ = "E",
        G = "G",
        P = "P",
        H = "H",
        C = "C",
        W = "W",
        M = "M",
        YFLIV = "Y"),
    group = list(
        DNSTA = c("D", "N", "S", "T", "A"), 
        EKRQ = c("E", "K", "R", "Q"),
        G = "G",
        P = "P",
        H = "H",
        C = "C",
        W = "W",
        M = "M",
        YFLIV = c("Y", "F", "L", "I", "V"))
) 

##  structure alignment
structure_alignments_Mirny <- list(
    color = c(
        DE = "#0082c8",
        KR = "#fabebe",
        NQST = '#aa6e28',
        GP = "#46f0f0",
        HWYF = "#f032e6",
        ACMLIV = "#e6194b"),
    symbol = c(
        DE = "D",
        KR = "K",
        NQST = "N",
        GP = "G",
        HWYF = "H",
        ACMLIV = "A"),
    group = list(
        DE = c("D","E"),
        KR = c("K","R"),
        NQST = c("N", "Q", "S", "T"),
        GP = c("G", "P"),
        HWYF = c("H", "W", "Y", "F"),
        ACMLIV = c("A", "C", "M", "L", "I", "V"))
) 

sequence_alignment_Dayhoff <- list(
    color = c(
        DENQ = "#0082c8",
        KRH = "#fabebe",
        STGPA = "#808000",
        C = '#3cb44b',
        WYF = "#808080",
        MLIV = "#e6beff"),
    symbol = c(
        DENQ = "D",
        KRH = "K",
        STGPA = "S",
        C = "C",
        WYF = "W",
        MLIV = "M"),
    group = list(
        DENQ = c("D","E", "N", "Q"),
        KRH = c("K","R", "H"),
        STGPA  = c("S", "T", "G", "P", "A"),
        C = "C",
        WYF = c("W", "Y", "F"),
        MLIV = c("M", "L", "I", "V"))
) 

chemistry_property_Mahler <- list(
    color = c(
        DE = "#0082c8",
        KRH = "#fabebe",
        QN = "#800000",
        ST = "#808000",
        P = "#00FF00",
        CM = '#3cb44b',
        WYF = "#808080",
        GALIV = "#46f0f0"),
    symbol = c(
        DE = "D",
        KRH = "K",
        QN = "Q",
        ST = "S",
        P = "P",
        CM = "C",
        WYF = "W",
        GALIV = "G"),
    group = list(
        DE= c("D","E"),
        KRH = c("K","R", "H"),
        QN = c("Q", "N"),
        ST  = c("S", "T"),
        P = "P",
        CM = c("C", "M"),
        WYF = c("W", "Y", "F"),
        GALIV = c("G", "A", "M", "L", "I", "V"))
)


cachedEnv$no <- no
cachedEnv$classic <- classic
cachedEnv$charge <- charge
cachedEnv$chemistry <- chemistry
cachedEnv$hydrophobicity <- hydrophobicity
cachedEnv$chemistry_property_Mahler<- chemistry_property_Mahler
cachedEnv$contact_potential_Maiorov <- contact_potential_Maiorov
cachedEnv$sequence_alignment_Dayhoff <- sequence_alignment_Dayhoff
cachedEnv$protein_blocks_Rogov <- protein_blocks_Rogov
cachedEnv$structure_alignments_Mirny <- structure_alignments_Mirny
cachedEnv$BLOSM50_L1 <- BLOSM50_L1
cachedEnv$BLOSM50_L2 <- BLOSM50_L2
cachedEnv$BLOSM50_L3 <- BLOSM50_L3
cachedEnv$BLOSM50_L4 <- BLOSM50_L4
cachedEnv$BLOSM50_L5 <- BLOSM50_L5
cachedEnv$BLOSM50_L6 <- BLOSM50_L6
cachedEnv$BLOSM50_L7 <- BLOSM50_L7
cachedEnv$BLOSM50_L8 <- BLOSM50_L8







