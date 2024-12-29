# fortran-variable-transformation

Fortran pure functions to perform Yeo-Johnson transformation [1] with MLE (maximum likelihood estimator) of the parameter ($\lambda$).

Golden Section Search [2] is used for fast ML estimation.

## Yeo-Johnson Transformation (YJ)

YJ is an extended version of Box-Cox transformation (YJ can handle negative values).

For a given value \( x \):

### If \( \lambda \neq 0 \):

\[
\phi(x; \lambda) =
\begin{cases} 
\frac{(x + 1)^\lambda - 1}{\lambda}, & \text{if } x \geq 0, \\
\frac{-((-x + 1)^{2 - \lambda} - 1)}{2 - \lambda}, & \text{if } x < 0.
\end{cases}
\]

### If \( \lambda = 0 \):

\[
\phi(x; \lambda) =
\begin{cases}
\ln(x + 1), & \text{if } x \geq 0, \\
-\ln(-x + 1), & \text{if } x < 0.
\end{cases}
\]

## References

1. Yeo, I. K., & Johnson, R. A. (2000). A new family of power
transformations to improve normality or symmetry. Biometrika.

2. Kiefer, J. (1953). Sequential minimization procedures.
Mathematical Statistics (pp. 145–172). Publisher: Wiley, New York.
