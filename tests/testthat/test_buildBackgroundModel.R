test_that("test_bg_p_bg", {
    bg <- buildBackgroundModel(seq.example, bg = "inputSet", 
        permutationSize = 2L)
    bg <- buildBackgroundModel(seq.example, bg = "nonInputSet", 
        proteome = proteome.example, permutationSize = 2L)
})

test_that("test_bg_p_combinations", {
    bg <- buildBackgroundModel(seq.example, bg = "inputSet", 
        model = "anchored", permutationSize = 2L)
})

test_that("test_bg_p_model", {
    bg <- buildBackgroundModel(seq.example, targetPosition = "Nterminus", 
        proteome = proteome.example, permutationSize = 2L)
    bg <- buildBackgroundModel(seq.example, targetPosition = "Cterminus", 
        proteome = proteome.example, permutationSize = 2L)
})

test_that("test_bg_p_others", {
    bg <- buildBackgroundModel(seq.example, uniqueSeq = FALSE, 
        proteome = proteome.example, permutationSize = 2L)
    bg <- buildBackgroundModel(seq.example, replacement = TRUE, 
        proteome = proteome.example, permutationSize = 2L)
})

