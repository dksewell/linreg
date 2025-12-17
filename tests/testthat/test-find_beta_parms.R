
test_that("Find beta parameters",{
  
  expect_no_error(
    beta_parms <-
      find_beta_parms(2/5,0.68,0.9)
  )
  expect_equal(
    beta_parms[1] / sum(beta_parms),
    c(shape1 = 2/5),
    tolerance = 0.05
  )
  expect_equal(
    pbeta(0.68,beta_parms[1],beta_parms[2]),
    0.9,
    tolerance = 0.05
  )
  
})