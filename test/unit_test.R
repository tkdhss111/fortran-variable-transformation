rm(list = ls())

#=====================================
# Test: Yeo-Johnson Tranformation
#
yj    <- read.csv("../build/test/yj.csv")
yj.l  <- read.csv("../build/test/yj_lambda.csv")
yj.lr <- read.csv("../build/test/yj_lambda_r.csv")

head(yj)
head(yj.l)
head(yj.lr)

identical(yj$transformed, yj.l$transformed)

# Compare results with R package: bestNormalize
if (!require(bestNormalize)) install.packages("bestNormalize")
yj.r <- yeojohnson(yj$log_normal, standardize = F)
yj.r
(diff <- sum(yj.r$x.t - yj$transformed))
# N.B. R-estimated lambda (-1.006137) is slightly different from mine (-1.00446498) due to the difference in optimization algorithm.

# yj.lr with lambda = -1.006137 (same as R)
(diff <- sum(yj.r$x.t - yj.lr$transformed))

# Histgram
hist(yj$std_normal,  breaks = 12)
hist(yj$log_normal,  breaks = 12)
hist(yj$transformed, breaks = 12)
hist(yj.r$x.t,       breaks = 12) # Same as the above.

#=====================================
# Test: Box-Cox Tranformation
#
bc    <- read.csv("../build/test/bc.csv")
bc.l  <- read.csv("../build/test/bc_lambda.csv")
bc.lr <- read.csv("../build/test/bc_lambda_r.csv")

head(bc)
head(bc.l)
head(bc.lr)

identical(bc$transformed, bc.l$transformed)

# Compare results with R package: bestNormalize
bc.r <- boxcox(bc.lr$log_normal, standardize = F)
bc.r
(diff <- sum(bc.r$x.t - bc$transformed))
# N.B. R-estimated lambda (-0.1727399) is slightly different from mine (-0.170742929) due to the difference in optimization algorithm.

# bc.lr with lambda = -0.1727399 (same as R)
(diff <- sum(bc.r$x.t - bc.lr$transformed))

hist(bc$std_normal,  breaks = 12)
hist(bc$log_normal,  breaks = 12)
hist(bc$transformed, breaks = 12)
hist(bc.r$x.t,       breaks = 12) # Histgram is as same as the above.

