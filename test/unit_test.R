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
# N.B. R-estimated lambda (-1.006137) is slightly different from mine (-1.00446498) due to the optimization algorithm.

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
#bc   <- read.csv("../build/test/bc.csv")
bc.lr <- read.csv("../build/test/bc_lambda_r.csv")
head(bc.lr)

bc.r <- boxcox(bc.lr$log_normal, standardize = F)
bc.r

# bc.lr with lambda = -0.1727399 (same as R)
(diff <- sum(yj.r$x.t - yj.lr$transformed))

hist(bc$std_normal,  breaks = 12)
hist(bc$log_normal,  breaks = 12)
hist(bc$transformed, breaks = 12)
hist(bc.r$x.t,       breaks = 12) # Histgram is as same as the above.

