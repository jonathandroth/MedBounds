
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MedBounds

<!-- badges: start -->
<!-- badges: end -->

The MedBounds package implements the methodology from the paper
[“Testing
Mechanisms”](https://www.jonathandroth.com/assets/files/TestingMechanisms_Draft.pdf)
by Soonwoo Kwon and Jonathan Roth. The package provides tests for the
“sharp null of full mediation”, which conjectures that the effect of a
treatment operates through a particular conjectured mechanism (or set of
mechanisms) M. It also provides lower bounds on the fraction of
“always-takers” who are affected by the treatment despite having the
same value of M regardless of treatment status.

## Installation

You can install the development version of MedBounds from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jonathandroth/MedBounds")
```

## Application to Baranov et al

We illustrate how the package can be used by walking through how the
code can be applied to the application of Baranov et al in Section XX of
the “Testing Mechanisms”. In Baranov, $D$ is a treatment for depression
and $Y$ is an index of outcomes for women’s financial empowerment. We
are interested in whether the effect of $D$ on $Y$ can be explained by a
mediator, or set of mediators, $M$. We consider three choices for $M$:
(a) the presence of a grandmother in the home, (b)relationship quality
with the woman’s husband, and (c) the combination of these two
mechanisms.

We first load the package and the data.

``` r
library(MedBounds)
## XX put data in package and load data here
```

We begin with the case where $M$ is a binary indicator for the presence
of a grandmother in the home.

XX test of sharp null

— XX Explain what significnce means — XX discuss discretization

XX lower bound on frac affected

XX density plot

By default, MedBounds imposes the monotonicity assumption that the
treatment can only increase the value of M. In this setting, this means
that everyone who would have a grandmother present without receiving CBT
treatment would also have one present when receiving CBT treatment. We
can relax this assumption by setting the XX parameter to be non-zero,
which bounds the number of “defiers” by XX.

XX redo the above with some interesting choice of defier bounds.

We next turn to the setting, where we are interested in testing whether
the effect is mediated by relationship quality with the husband, which
is measured on a 1-5 scale. We can again test the sharp null and
estimate a lower bound on the fraction affected.

XX XX

Finally, we can test the null hypothesis that the treatment effect is
explained by the combination of the two mechanisms.
