## Supplementary Material

### Manuscript

The current repository provides the Supplementary Material for the
manuscript

> Zietkiewicz P and Kosmidis I (2023). Bounded-memory adjusted scores
> estimation in generalized linear models with large data sets.

### Contents

The Supplementary Material provides

i) the Supplementary Material document that is cross-referenced in the
manuscript and contains all numerical results and figures from the
case study of Section 4.

ii) R code to reproduce all numerical results in the main text and in
the Supplementary Materials document.

The code is organized in the two directories `diverted-flights` and
`biglm`.

The `diverted-flights` directory, apart form code and data, provides a
`README` file with specific instructions to reproduce the numerical
results, along with the specific versions of the contributed R
packages that have been used.

The `biglm` directory has a port of the [**biglm** R
package](https://cran.r-project.org/package=biglm), which implements
the one- and two-pass IWLS variants for solving the bias-reducing
adjusted score equations ([Firth,
1993](https://doi.org/10.1093/biomet/80.1.27)) and for maximum
Jeffreys'-penalized likelihood estimation ([Kosmidis & Firth,
2021](https://doi.org/10.1093/biomet/asaa052)).

### Authors

Patrick Zietkiewicz <patrick.zietkiewicz@warwick.ac.uk>

Ioannis Kosmidis <ioannis.kosmidis@warwick.ac.uk>

