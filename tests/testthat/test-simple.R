context("simple models as special cases of the general model")

test_that("two_parts", {
    # log-lin simple hurdle uncorrelated / correlated
    L100I <- mhurdle(shows ~ educ | linc, Interview, dist = "ln", h2 = FALSE, corr = FALSE,
                     check_gradient =TRUE, scaled = FALSE, robust = FALSE, start = c(rep(0.1, 4), 1))
    L100D <- update(L100I, corr = TRUE, start = c(rep(0.1, 4), 1, 0))
    expect_equal(L100I$value, L100D$value)
    expect_true(all(abs(L100D$gradient$rel_diff < 1E-05)))
    expect_equal(L100I$gradient, L100D$gradient[- 6, ])

    # log-lin simple / double hurdle uncorrelated
    # (note that the numerical derivative for mu is irrelevant for mu = 0)
    L110I <- update(L100I, h2 = TRUE, start = c(rep(0.1, 4), 1, 0))
    expect_equal(L100I$value, L100D$value)
    expect_equal(L100I$gradient, L110I$gradient[- 6, ])
    expect_true(all(abs(L110I$gradient$rel_diff[- 6]) < 1E-05))

    # lin simple hurdle uncorrelated / correlated
    N100I <- mhurdle(shows ~ educ | linc, Interview, dist = "n", h2 = FALSE, corr = FALSE,
                     check_gradient =TRUE, scaled = FALSE, robust = FALSE, start = c(rep(0.1, 4), 1))
    N100D <- update(N100I, corr = TRUE, start = c(rep(0.1, 4), 1, 0))
    expect_equal(N100I$value, N100D$value)
    expect_true(all(abs(N100D$gradient$rel_diff < 1E-05)))
    expect_equal(L100I$gradient, L100D$gradient[- 6, ])
})

test_that("one_equation", {
    K <- c(n = 3, ln = 4, bc = 5)
    for (dist in c("n", "ln")){
        one_eq <- mhurdle(shows ~ 0 | linc, Interview, dist = dist, h2 = TRUE, check_gradient = TRUE)
        expect_true(all(one_eq$gradient[, "rel_diff"] < 1E-05))
        one_eq <- mhurdle(shows ~ 0 | linc, Interview, dist = dist, h2 = TRUE, check_gradient = TRUE, start = rep(1, K[dist]))
        expect_true(all(one_eq$gradient[, "rel_diff"] < 1E-05))
    }
})
