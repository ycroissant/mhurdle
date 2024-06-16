# mhurdle 1.3-1

- For robust estimations (the default), the model is now internally
  updated with `robust = FALSE` and `iterlim = 0` so that the gradient
  and the hessian of the structural model are computed.
- Interface to `sandwich` (`estfun` and `bread` methods) and to `nonnest2`
  (`llcont` method).

# mhurdle 1.3-0

- `predict(object, newdata = data)` where `data` is the data frame used
  to fit the model now returns the same as `fitted(object)`. Bug
  fixed thanks to Achim Zeileis and Rebekka Topp.
- Interface with foreign packages (`prediction`, `margins`, and
  `modelsummary`).
- Unit tests added.

# mhurdle 1.1-6

- The EV formula is fixed for the log-normal model.
- Improved version of the `texreg` method.

# mhurdle 1.0-1

- Major revision for the code and the vignette; bc and ihs
  transformation, heteroscedasticity are introduced.

# mhurdle 0.1-3

- Minor changes of the vignette.

# mhurdle 0.1-2

- Major update of the whole code.
- Much improved vignette.

# mhurdle 0.1-0

- Some encoding problems in the vignette are fixed.
	
# mhurdle 0.1-0

- First CRAN release.
