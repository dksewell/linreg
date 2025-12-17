
test_that("Test lm_b method for getting posterior draws",{
  
  # Create data and fit lm_b model
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  expect_no_error(
    fita <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "conj")
  )
  expect_true(
    is.matrix(get_posterior_draws(fita,
                                  n_draws = 10))
  )
  
  
})