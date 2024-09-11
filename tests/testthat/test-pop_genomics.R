test_that("Population genomics functions work", {
  vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
  index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")

  expect_true(class(FixedSites(vcf_file)) == "integer")
  expect_true(class(FixedSites(vcf_file, window_size = 10000, skip_size = 5000)) == "data.frame")

  expect_true(class(SegregatingSites(vcf_file)) == "integer")
  expect_true(class(SegregatingSites(vcf_file, window_size = 10000, skip_size = 5000)) == "data.frame")

  expect_true(class(SingletonSites(vcf_file)) == "integer")
  expect_true(class(SingletonSites(vcf_file, window_size = 10000, skip_size = 5000)) == "data.frame")


  pop1_individuals <- c("tsk_0", "tsk_1", "tsk_2")
  pop2_individuals <- c("tsk_3", "tsk_4", "tsk_5")
  expect_true(class(PrivateAlleles(vcf_file, pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals)) == "numeric")
  expect_true(class(PrivateAlleles(vcf_file, pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals, window_size = 10000, skip_size = 5000)) == "data.frame")

  expect_true(class(Heterozygosity(vcf_file)) == "numeric")
  expect_true(class(Heterozygosity(vcf_file, window_size = 10000, skip_size = 5000)) == "data.frame")

  total_sequence_length <- 999299
  expect_true(class(Pi(vcf_file, seq_length = total_sequence_length)) == "numeric")
  expect_true(class(Pi(vcf_file, seq_length = total_sequence_length, window_size = 10000, skip_size = 5000)) == "data.frame")

  expect_true(class(TajimasD(vcf_file, seq_length = total_sequence_length)) == "numeric")
  expect_true(class(TajimasD(vcf_file, seq_length = total_sequence_length, window_size = 10000, skip_size = 5000)) == "data.frame")

  expect_true(class(WattersonsTheta(vcf_file, seq_length = total_sequence_length)) == "numeric")
  expect_true(class(WattersonsTheta(vcf_file, seq_length = total_sequence_length, window_size = 10000, skip_size = 5000)) == "data.frame")

  expect_true(class(Dxy(vcf_file, pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals, seq_length = total_sequence_length)) == "numeric")
  expect_true(class(Dxy(vcf_file, pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals, seq_length = total_sequence_length, window_size = 10000, skip_size = 5000)) == "data.frame")

  expect_true(class(Fst(vcf_file, pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals)) == "numeric")
  expect_true(class(Fst(vcf_file, pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals, window_size = 10000, skip_size = 5000)) == "data.frame")

  expect_true(class(OneDimSFS(vcf_file)) == "numeric")
  expect_true(class(TwoDimSFS(vcf_file, pop1_individuals = pop1_individuals, pop2_individuals = pop2_individuals))[1] == "matrix")

})
