## create a cached environment to store those package-wide accessible constants.
cacheEnv <- new.env(parent = emptyenv())

## color coding for amino acids
auto <- c('A' = '#CCFF00',
          'C' = '#FFFF00',
          'D' = '#FF0000',
          'E' = '#FF0066',
          'F' = '#00FF66',
          'G' = '#FF9900',
          'H' = '#0066FF',
          'I' = '#66FF00',
          'K' = '#6600FF',
          'L' = '#33FF00',
          'M' = '#00FF00',
          'N' = '#CC00FF',
          'P' = '#FFCC00',
          'Q' = '#FF00CC',
          'R' = '#0000FF',
          'S' = '#FF3300',
          'T' = '#FF6600',
          'V' = '#99FF00',
          'W' = '#00CCFF',
          'Y' = '#00FFCC')

namehash <- c(
    "Ala" = "A",
    "Arg" = "R",
    "Asn" = "N",
    "Asp" = "D",
    "Cys" = "C",
    "Glu" = "E",
    "Gln" = "Q",
    "Gly" = "G",
    "His" = "H",
    "Ile" = "I",
    "Leu" = "L",
    "Lys" = "K",
    "Met" = "M",
    "Phe" = "F",
    "Pro" = "P",
    "Ser" = "S",
    "Thr" = "T",
    "Trp" = "W",
    "Tyr" = "Y",
    "Val" = "V")

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
    group = list(
        nonpolar_aliphatic = c("A", "G", "L", "M", "I", "V"),
        polar_uncharged = c("C", "P", "Q", "S", "T"),
        aromatic = c("F", "W", "Y"),
        positively_charged = c("H", "K", "N", "R"),
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
        hydrophobic = c("A", "F", "I", "L", "M", "P", "V", "W"),
        polar = c("C", "G", "S", "T", "Y"),
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
        hydrophilic = c("D", "E", "K", "N", "Q", "R"),
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

cacheEnv$auto <- auto
cacheEnv$namehash <- namehash
cacheEnv$classic <- classic
cacheEnv$chemistry <- chemistry
cacheEnv$hydrophobicity <- hydrophobicity
cacheEnv$charge <- charge
