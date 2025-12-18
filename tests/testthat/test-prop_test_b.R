

test_that("Test prop_test_b for a single population",{
  
  # Test input methods
  expect_no_error(
    fita <-
      prop_test_b(14,
                  19)
  )
  expect_no_error(
    fitb <-
      prop_test_b(14,
                  n_total = 14 + 19)
  )
  expect_equal(fita[c(1:5,7)],
               fitb[c(1:5,7)])
  
  # Test output
  expect_equal(fita$successes,14)
  expect_equal(fita$failures,19)
  expect_equal(fitb$failures,19)
  expect_type(fita$posterior_mean,"double")
  expect_type(fita$CI,"double")
  expect_type(fita$PI,"double")
  expect_type(fita$posterior_parameters,"double")
  expect_s3_class(fita$prop_plot,
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # Test probability of comparison
  expect_no_error(
    fitc <-
      prop_test_b(14,
                  19,
                  p = 0.45)
  )
  expect_type(fitc$Pr_less_than_p,"double")
  
})



test_that("Test prop_test_b for two populations",{
  
  # Test input methods
  expect_no_error(
    fita <-
      prop_test_b(c(14,22),
                  c(19,45))
  )
  expect_no_error(
    fitb <-
      prop_test_b(c(14,22),
                  n_total = c(14,22) + c(19,45))
  )
  expect_equal(fita[c(1:10,12)],
               fitb[c(1:10,12)])
  
  # Test output
  expect_equal(fita$successes,c(14,22))
  expect_equal(fita$failures,c(19,45))
  expect_equal(fitb$failures,c(19,45))
  expect_type(fita$posterior_mean_pop1,"double")
  expect_type(fita$posterior_mean_pop2,"double")
  expect_type(fita$CI_pop1,"double")
  expect_type(fita$CI_pop2,"double")
  expect_type(fita$PI_pop1,"double")
  expect_type(fita$PI_pop2,"double")
  expect_type(fita$posterior_parameters$population_1,"double")
  expect_type(fita$posterior_parameters$population_2,"double")
  expect_s3_class(fita$prop_plot,
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # Test probability of comparison doesn't do anything
  expect_no_error(
    fitc <-
      prop_test_b(c(14,22),
                  c(19,45),
                  p = 0.45)
  )
  expect_equal(fita[c(1:10,12)],
               fitc[c(1:10,12)])
  
})