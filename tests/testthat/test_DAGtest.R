test_that("test_testDAU_group", {
    bg <- buildBackgroundModel(seq.example, proteome = proteome.example, 
        permutationSize = 10L)
    t <- testDAU(seq.example, bg)
    expect_true(t@zscore["K", "AA0"] > 100)
    expect_equal(0, t@pvalue["K", "AA0"])
    t <- testDAU(seq.example, bg, group = "classic")
    expect_true(t@zscore["positively_charged", "AA0"] > 50)
    expect_equal(0, t@pvalue["positively_charged", "AA0"])
})

