
test_that("Test semi-parametric survival curve fitting",{
  set.seed(2025)
  N = 300
  test_data = 
    data.frame(outcome = 
                 rweibull(N,2,5))
  test_data$observed = 
    ifelse(test_data$outcome >= 7, 0, 1)
  test_data$outcome =
    ifelse(dplyr::near(test_data$observed,1), test_data$outcome, 7)
  
  # Check fitting function
  expect_no_error(
    fita <-
      survfit_b(Surv(test_data$outcome,
                     test_data$observed))
  )
  
  # Check print
  expect_no_error(
    fita
  )
  
  # Check plotting
  expect_no_error(
    survplot <- 
      plot(fita)
  )
  expect_type(survplot,"list")
  expect_s3_class(survplot$plot,
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(survplot$data,c("tbl_df", "tbl", "data.frame"))
  
  
})
