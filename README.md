## Overview

Data collected by means of households' expenditure survey may present
a large proportion of zero expenditures due to many households
recording, for one reason or another, no expenditure for some
items. Since the seminal paper of Tobin (1958), a large econometric
literature has been developed to deal correctly with this problem of
zero observations. In particular, a good selection mechanism was
introduced by Cragg (1971) and a purchasing mechanism by Deaton and
Irish (1984). We propose an encompassing approach with a general three
equations model for which a zero expense can be observed either
because:

- the good is not selected,
- the good is selected but not consumed because of lack of financial
  resources,
- the good is selected and consumed, but purchased infrequently so
  that no expense is observed during the period of the survey.

`mhurdle` provides a set of tools to estimate (by maximum likelihood)
 and test (using especially vuong test) this generalized hurdle model.

## Installation

`mhurdle` is  on `CRAN`.

```
install.packages("mhurdle")
```

For the development version, use 

```
install.packages("devtools")
devtools::install_github("ycroissant/mhurdle")
```
