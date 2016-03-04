require("dagLogo") || stop("unable to load Package:dagLogo")
require("Biostrings") || stop("unable to load Package:Biostrings")
data(seq.example)
data(proteome.example)
library(testthat)

test_check("dagLogo")