Changes since version 1.3-0
  * for robust estimations (the default), the model is now internally
    updated with robust = FALSE and iterlim = 0 so that the gradient
    and the hessian of the structural model are computed 
  * interface to sandwich (estfun and bread methods) and to nonnest2
    (llcont method)

Changes since version 1.2-0

  * predict(object, newdata = data) where data is the data.frame used
      to fit the model now returns the same as fitted(object) bug
      fixed thanks to Achim Zeileis and Rebekka Topp
  * interface with foreign packages (prediction, margins and
      modelsummary)
  * unit tests added

Changes since version 1.1-7

  * the EV formula is fixed for the log-normal model
  * improved version of the texreg method

Changes since version 0.1-4

  * major revision for the code and the vignette ; bc and ihs
    transformation, heteroscedasticity are introduced

Changes since version 0.1-2 :

  * minor changes of the vignette

Changes since version 0.1-1 :

  * major update of the whole code
  * much improved vignette
	
Changes since version 0.1-0 :
	
  * some encoding problems in the vignette are fixed
	
