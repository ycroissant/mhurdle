context("analytical vs numerical gradient / hessian")


test_that("gradients", {
    data("Interview", package = "mhurdle")
    for (form in list(shows ~ educ | linc | age, shows ~ educ | linc, shows ~ 0 | linc | age)){
        for (dist in c("n", "ln")){
            for (h2 in c(TRUE, FALSE)){
                for (corr in c(TRUE, FALSE)){
                    za <- mhurdle(form, data = Interview, h2 = h2, corr = corr, dist = dist, check_gradient = TRUE)$gradient
                    cat(paste(deparse(form), ", dist = ", dist, ", h2 = ", h2, ", corr = ", corr, "\n"))
                    expect_true(all(za[, "rel_diff"] < 1E-04))
                }
            }
        }
    }
})


