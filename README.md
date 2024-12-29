# fortran-variable-transformation (work-in-progress)

Fortran pure functions to perform Yeo-Johnson transformation [1] with MLE (maximum likelihood estimator) of the parameter ($\lambda$).

Golden Section Search [2] is used for fast ML estimation.

## References

1. Yeo, I. K., & Johnson, R. A. (2000). A new family of power
transformations to improve normality or symmetry. Biometrika.

2. Kiefer, J. (1953). Sequential minimization procedures.
Mathematical Statistics (pp. 145–172). Publisher: Wiley, New York.

# Yeo-Johnson Transformation [by ChatGPT]

The **Yeo-Johnson transformation** is a family of power transformations used to normalize data. It can handle both positive and negative values, unlike the Box-Cox transformation, which only works for strictly positive values.

## Formula

For a given value \( x \), the Yeo-Johnson transformation is defined as:

### If \( \lambda \neq 0 \):
\[
y(x; \lambda) =
\begin{cases} 
\frac{(x + 1)^\lambda - 1}{\lambda}, & \text{if } x \geq 0, \\
\frac{-((-x + 1)^{2 - \lambda} - 1)}{2 - \lambda}, & \text{if } x < 0.
\end{cases}
\]

### If \( \lambda = 0 \):
\[
y(x; \lambda) =
\begin{cases}
\ln(x + 1), & \text{if } x \geq 0, \\
-\ln(-x + 1), & \text{if } x < 0.
\end{cases}
\]

## Explanation

- \( x \) is the input value, and \( \lambda \) is the transformation parameter.
- When \( \lambda = 0 \), the transformation reduces to the logarithmic transformation.
- The transformation is applied differently for positive and negative values of \( x \):
  - For \( x \geq 0 \), the transformation is applied using the first formula.
  - For \( x < 0 \), the transformation uses a different formula that ensures the output remains well-defined.

## Purpose

The Yeo-Johnson transformation is useful in making data more normally distributed, which is often a prerequisite for many statistical methods and machine learning algorithms. It is particularly useful when your data contains both positive and negative values.
