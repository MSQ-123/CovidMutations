context("test basic output figure of AssayMutRatio")



test_that("test basic figure",{
  td<-tempdir()
  #setTmpDir(td)

  Total <- 11000
  data("nucmerr")
  data("assays")
  AssayMutRatio(nucmerr = nucmerr,
                assays = assays,
                totalsample = Total,
                plotType = "logtrans",
                outdir = td)

  expect_true(file.exists(file.path(td,"log2_Mutation_Ratio.png")))



})
