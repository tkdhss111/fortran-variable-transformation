rm(list = ls())

d1 <- read.csv("../build/test/data.csv")
d2 <- read.csv("../build/test/data_user_lambda.csv")

head(d1)
head(d2)

# Test: User-defined lambda
identical(d1$yj_transformed, d2$yj_transformed)

# Compare results with R package: bestNormalize
if (!require(bestNormalize)) install.packages("bestNormalize")
z.yj.r <- yeojohnson(d1$log_normal, standardize = F)
z.yj.r
# N.B. R-estimated lambda (-1.006137) is slightly different from mine (-1.00446498) due to the optimization algorithm.
(diff <- sum(z.yj.r$x.t - d1$yj_transformed))

hist(d1$std_normal,     breaks = 12)
hist(d1$log_normal,     breaks = 12)
hist(d1$yj_transformed, breaks = 12)
hist(z.yj.r$x.t,        breaks = 12) # Histgram is as same as the above.
