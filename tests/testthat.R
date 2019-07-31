library(testthat)
library(OOSAP)

test_check("OOSAP", reporter = "progress")

# if (requireNamespace("lintr", quietly = TRUE)) {
#     context("lints")
#     test_that("Package Style", {
#         lintr::expect_lint_free()
#     })
# }