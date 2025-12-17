
test_that("find inverse gamma parameters",{
  
  # Test method aimed at linear regression
  hypothetical_s2_y = 2.0
  lower_R2 = 0.05
  upper_R2 = 0.85
  
  expect_no_error(
    igamma_parms <- 
      find_invgamma_parms(response_variance = hypothetical_s2_y,
                          lower_R2 = lower_R2,
                          upper_R2 = upper_R2,
                          probability = 0.8)
  )
  
  expect_equal(
    extraDistr::pinvgamma(hypothetical_s2_y * (1.0 - upper_R2),
                          igamma_parms[1] * 0.5,
                          igamma_parms[2] * 0.5),
    0.1,
    tolerance = 0.025)
  expect_equal(
    extraDistr::pinvgamma(hypothetical_s2_y * (1.0 - lower_R2),
                          igamma_parms[1] * 0.5,
                          igamma_parms[2] * 0.5,
                          lower.tail = FALSE),
    0.1,
    tolerance = 0.025)
  
  
  # Test more arbitrary method
  igamma_parms2 <-
    find_invgamma_parms(hypothetical_s2_y * (1.0 - upper_R2),
                        hypothetical_s2_y * (1.0 - lower_R2),
                        probability = 0.8)
  expect_equal(igamma_parms,
               igamma_parms2)
  
})