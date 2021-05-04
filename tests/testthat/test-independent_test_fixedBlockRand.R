context("Independent testing for fixedBlockRand.R")

test_that("fixedBlockRand returns an appropriate object: a vector of length n of treatment group",{
  s <- fixedBlockRand(n=10)
  expect_equal(length(s),10)
})

test_that("the final block is a partial block if `n` is not an exact multiple of the block length",{
  n <- 10
  block <- c(0,0,1,1)
  s <- fixedBlockRand(n=n, block=block)
  lastblock <- s[(length(block)*floor(n/length(block))+1):length(s)]
  
  freq <- table(block)
  freq_lastblock <- table(lastblock)
  
  if (0 %in% lastblock) expect_lte(as.numeric(freq_lastblock[names(freq_lastblock)==0]),as.numeric(freq[names(freq)==0]))
  if (1 %in% lastblock) expect_lte(as.numeric(freq_lastblock[names(freq_lastblock)==1]),as.numeric(freq[names(freq)==1]))
})
