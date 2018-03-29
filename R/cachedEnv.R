## create a cached environment to store those package-wide accessible constants.
cachedEnv <- new.env(parent = emptyenv())

## color and character symbol encodings and grouping for amino acids
no <- list(
    color = c('A' = '#CCFF00', 'C' = '#FFFF00', 'D' = '#FF0000',
              'E' = '#FF0066', 'F' = '#00FF66', 'G' = '#FF9900',
              'H' = '#0066FF', 'I' = '#66FF00', 'K' = '#6600FF',
              'L' = '#33FF00', 'M' = '#00FF00', 'N' = '#CC00FF',
              'P' = '#FFCC00', 'Q' = '#FF00CC', 'R' = '#0000FF',
              'S' = '#FF3300', 'T' = '#FF6600', 'V' = '#99FF00',
              'W' = '#00CCFF', 'Y' = '#00FFCC'),
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
        negatively_charged = "#FFB32C"
    ),
    symbol = c(
        nonpolar_aliphatic = "M",
        polar_uncharged = "U",
        aromatic = "A",
        positively_charged = "P",
        negatively_charged = "N"
    ),
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
        negatively_charged = c("D", "E")
    )
)

chemistry <- list(
    color = c(
        hydrophobic = "#000000",
        polar = "#00811B",
        basic = "#2000C7",
        neutral  = "#800080",
        acidic = "#D00001"
    ),
    symbol = c(
        hydrophobic = "H",
        polar = "P",
        basic = "B",
        neutral  = "U",
        acidic = "A"
    ),
    group = list(
        hydrophobic = c("G", "A", "I", "L", "V", "M", "P", "F","W"),
        polar = c("C", "S", "T", "Y"),
        basic = c("H", "K", "R"),
        neutral = c("N", "Q"),
        acidic = c("D", "E")
    )
)

hydrophobicity <- list(
    color = c(
        hydrophilic = '#000000',
        neutral = '#00811B',
        hydrophobic  = '#2000C7'
    ),
    symbol = c(
        hydrophilic = 'I',
        neutral = 'U',
        hydrophobic = 'O'
    ),
    group = list(
        hydrophilic = c("D", "E", "K", "R", "N", "Q"),
        neutral = c("A", "G", "H", "P", "S", "T"),
        hydrophobic = c("C", "F", "I", "L", "M", "V", "W", "Y")
    )
)

charge <- list(
    color = c(
        positive = "#FFB32C",
        neutral = "#2000C7",
        negative = "#CCCCCC"
    ),
    symbol = c(
        positive = "P",
        neutral = "U",
        negative = "N"
    ),
    group = list(
        positive = c("H", "K", "R"),
        neutral = c("A", "C", "F", "G", "I", "L", "M", "N", "P", 
                    "Q", "S", "T", "V","W","Y"),
        negative = c("D", "E")
    )
) 

## 8-level schemes for amino acid alphabet reduction derived from correlations 
## based on the BLOSUM50 similarity matrix. More schemes are available from
## Stephenson and Freeland. Unwearthing the root of amino acid similarity. 
## J Mol. Evol. 2013 (77):159-169. Smith and Smith. Protein Eng. 1992(5):35-41.
## And many other literature.

BLOSM50_L1 <- list(
    color = c(
        LVIMCAGSTPFYW = "#33FF00",
        EDNQKRH = "#FF0066"
    ),
    symbol = c(
        LVIMCAGSTPFYW = "L",
        EDNQKRH = "E"
    ),
    group = list(
        LVIMCAGSTPFYW = c("L", "V", "I", "M", "C", "A", "G", "S", "T", "P", 
                          "F", "Y", "W"),
        EDNQKRH = c("E", "D", "N", "Q", "K", "R", "H")
    )
) 

BLOSM50_L2 <- list(
    color = c(
        LVIMCAGSTP = "#33FF00",
        FYW = "#00FF66",
        EDNQKRH = "#FF0066"
    ),
    symbol = c(
        LVIMCAGST = "L",
        FYW = "F",
        EDNQKRH = "E"
    ),
    group = list(
        LVIMCAGSTP = c("L", "V", "I", "M", "C", "A", "G", "S", "T", "P"), 
        FYW = c("F", "Y", "W"),
        EDNQKRH = c("E", "D", "N", "Q", "K", "R", "H")
    )
) 

BLOSM50_L3 <- list(
    color = c(
        LVIMC = "#33FF00",
        AGSTP = "#CCFF00",
        FYW = '#00FF66',
        EDNQKRH = "#FF0066"
    ),
    symbol = c(
        LVIMC = "L",
        AGSTP = "A",
        FYW = "F",
        EDNQKRH = "E"
    ),
    group = list(
        LVIMC = c("L", "V", "I", "M", "C"), 
        AGSTP = c("A", "G", "S", "T", "P"),
        FYW = c("F", "Y", "W"),
        EDNQKRH = c("E", "D", "N", "Q", "K", "R", "H")
    )
) 

BLOSM50_L4 <- list(
    color = c(
        LVIMC = "#33FF00",
        AGSTP = "#CCFF00",
        FYW = '#00FF66',
        EDNQ = "#FF0066",
        KRH = "#6600FF"
    ),
    symbol = c(
        LVIMC = "L",
        AGSTP = "A",
        FYW = "F",
        EDNQ = "E",
        KRH = "K"
    ),
    group = list(
        LVIMC = c("L", "V", "I", "M", "C"), 
        AGSTP = c("A", "G", "S", "T", "P"),
        FYW = c("F", "Y", "W"),
        EDNQ = c("E", "D", "N", "Q"),
        KRH = c("K", "R", "H")
    )
) 

BLOSM50_L5 <- list(
    color = c(
        LVIMC = "#33FF00",
        AGST = "#CCFF00",
        P = "#FFCC00",
        FYW = "#00FF66",
        EDNQ = "#FF0066",
        KRH = "#6600FF"
    ),
    symbol = c(
        LVIMC = "L",
        AGSTP = "A",
        P = "P",
        FYW = "F",
        EDNQ = "E",
        KRH = "K"
    ),
    group = list(
        LVIMC = c("L", "V", "I", "M", "C"), 
        AGST = c("A", "G", "S", "T"),
        P = "P",
        FYW = c("F", "Y", "W"),
        EDNQ = c("E", "D", "N", "Q"),
        KRH = c("K", "R", "H")
    )
) 


BLOSM50_L6 <- list(
    color = c(
        LVIMC = "#33FF00",
        AG = "#CCFF00",
        ST = "#FF3300",
        P = "#FFCC00",
        FYW = "#00FF66",
        EDNQ = "#FF0066",
        KR = "#6600FF",
        H = "#0066FF"
    ),
    symbol = c(
        LVIMC = "L",
        AG = "A",
        ST = "S",
        P = "P",
        FYW = "F",
        EDNQ = "E",
        KR = "K",
        H = "H"
    ),
    group = list(
        LVIMC = c("L", "V", "I", "M", "C"), 
        AG = c("A", "G"),
        ST = c("S", "T"),
        P = "P",
        FYW = c("F", "Y", "W"),
        EDNQ = c("E", "D", "N", "Q"),
        KR = c("K", "R"),
        H = "H"
    )
) 

BLOSM50_L7 <- list(
    color = c(
        LVIM = "#33FF00",
        C = "#FFFF00",
        A = "#CCFF00",
        G = "#FF9900",
        ST = "#FF3300",
        P = "#FFCC00",
        FYW = "#00FF66",
        EDNQ = "#FF0066",
        KR = "#6600FF",
        H = "#0066FF"
    ),
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
        H = "H"
    ),
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
        H = "H"
    )
) 

BLOSM50_L8 <- list(
    color = c(
        LVIM = "#33FF00",
        C = "#FFFF00",
        A = "#2000C7",
        G = "#FF9900",
        S = "#FF3300",
        T = "#FF6600",
        P = "#FFCC00",
        FY = "#00FF66",
        W = "#00CCFF",
        E = "#FF0066",
        D = "#FF0000",
        N = "#CC00FF",
        Q = "#FF00CC",
        KR = "#6600FF",
        H = "#0066FF"
    ),
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
        H = "H"
    ),
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
        H = "H"
    )
) 


## Spatial frequency
contact_potential_Maiorov <- list(
    color = c(
        DENQ = "#FF0000",
        KR = "#6600FF",
        G = '#FF9900',
        P = "#FFCC00",
        AV = "#CCFF00",
        STHWY = "#FF3300",
        CFMLI = "#FFFF00"
    ),
    symbol = c(
        DENQ = "D",
        KR = "K",
        G = "G",
        P = "P",
        AV = "A",
        STHWY = "S",
        CFMLI = "C"
    ),
    group = list(
        DENQ = c("D", "E", "N", "Q"), 
        KR = c("K", "R"),
        G = "G",
        P = "P",
        AV = c("A", "V"),
        STHWY = c("S", "T", "H", "W", "Y"),
        CFMLI = c("C", "F", "M", "L", "I")
    )
) 

## protein blocks
protein_blocks_Rogov <- list(
    color = c(
        DNSTA = "#FF0000",
        EKRQ = "#FF0066",
        G = '#FF9900',
        P = "#FFCC00",
        H = "#0066FF",
        C = "#FFFF00",
        W = "#00CCFF",
        M = "#00FF00",
        YFLIV = "#00FFCC"
    ),
    symbol = c(
        DNSTA = "D",
        EKRQ = "E",
        G = "G",
        P = "P",
        H = "H",
        C = "C",
        W = "W",
        M = "M",
        YFLIV = "Y"
    ),
    group = list(
        DNSTA = c("D", "N", "S", "T", "A"), 
        EKRQ = c("E", "K", "R", "Q"),
        G = "G",
        P = "P",
        H = "H",
        C = "C",
        W = "W",
        M = "M",
        YFLIV = c("Y", "F", "L", "I", "V")
    )
) 

##  structure alignment
structure_alignments_Mirny <- list(
    color = c(
        DE = "#FF0000",
        KR = "#6600FF",
        NQST = '#CC00FF',
        GP = "#FF9900",
        HWYF = "#0066FF",
        ACMLIV = "#CCFF00"
    ),
    symbol = c(
        DE = "D",
        KR = "K",
        NQST = "N",
        GP = "G",
        HWYF = "H",
        ACMLIV = "A"
    ),
    group = list(
        DE = c("D","E"),
        KR = c("K","R"),
        NQST = c("N", "Q", "S", "T"),
        GP = c("G", "P"),
        HWYF = c("H", "W", "Y", "F"),
        ACMLIV = c("A", "C", "M", "L", "I", "V")
    )
) 

sequence_alignment_Dayhoff <- list(
    color = c(
        DENQ = "#FF0000",
        KRH = "#6600FF",
        STGPA = "#FF3300",
        C = '#FFFF00',
        WYF = "#00CCFF",
        MLIV = "#00FF00"
    ),
    symbol = c(
        DENQ = "D",
        KRH = "K",
        STGPA = "S",
        C = "C",
        WYF = "W",
        MLIV = "M"
    ),
    group = list(
        DENQ = c("D","E", "N", "Q"),
        KRH = c("K","R", "H"),
        STGPA  = c("S", "T", "G", "P", "A"),
        C = "C",
        WYF = c("W", "Y", "F"),
        MLIV = c("M", "L", "I", "V")
    )
) 

chemistry_property_Mahler <- list(
    color = c(
        DE = "#FF0000",
        KRH = "#6600FF",
        QN = "#FF00CC",
        ST = "#FF3300",
        P = "#FFCC00",
        CM = '#FFFF00',
        WYF = "#00CCFF",
        GALIV = "#FF9900"
    ),
    symbol = c(
        DE = "D",
        KRH = "K",
        QN = "Q",
        ST = "S",
        P = "P",
        CM = "C",
        WYF = "W",
        GALIV = "G"
    ),
    group = list(
        DE= c("D","E"),
        KRH = c("K","R", "H"),
        QN = c("Q", "N"),
        ST  = c("S", "T"),
        P = "P",
        CM = c("C", "M"),
        WYF = c("W", "Y", "F"),
        GALIV = c("G", "A", "M", "L", "I", "V")
    )
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







