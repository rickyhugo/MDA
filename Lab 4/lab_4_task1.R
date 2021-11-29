## ---- warning=FALSE-------------------------------------------------------------------------------------
# import libraries
library(readxl)
library(stats)
library(matlib)
library(ggplot2)
library(qqplotr)
library(broom)
library(Hmisc)
library(dplyr)


## -------------------------------------------------------------------------------------------------------
# read data
d = as.data.frame(read_excel('Knee.xlsx'))
d[d == "M"] = NA
# summary(d)


## -------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Task 1

# Variable definitions:
# I_Q_conNm_weight (injured, concentric)
# C_Q_conNm_weight (healthy, concentric) 
# I_Q_eccNm_weight (injured, eccentric)
# C_Q_eccNm_weight (healthy, eccentric)

# Test: Multivariate Hotelling's T-squared test statistic (Paired comparisons)

# extract relevant variables
d_t1 = cbind(
  group=d$Grupp,
  I_Q_con=d$I_Q_conNm_weight,
  I_Q_ecc=d$I_Q_eccNm_weight,
  I_H_con=d$I,
  I_H_ecc=d$I_H_eccNm_weight,
  C_Q_con=d$C_Q_conNm_weight,
  C_Q_ecc=d$C_Q_eccNm_weight,
  C_H_con=d$C_H_conNm_weight,
  C_H_ecc=d$C_H_eccNm_weight
)

# remove NA's and convert extracted data to data frame
D_t1_chr = as.data.frame(na.omit(d_t1))
D_t1 = as.data.frame(sapply(D_t1_chr, as.numeric))

# split into separate data frame for each group
D_t1_gr1 = D_t1[D_t1$group == 1, ]
D_t1_gr2 = D_t1[D_t1$group == 2, ]


## -------------------------------------------------------------------------------------------------------
# compute paired comparison statistics
paired_comp = function(gr1, gr2, a=.05) {
  # compute parameters
  D = gr1 - gr2
  D_bar = as.matrix(colMeans(D))
  S = cov(D)
  n = dim(gr1)[1]
  p = dim(gr1)[2]
  
  # compute statistics (H0: delta = 0)
  H_T2 = n * t(D_bar) %*% inv(S) %*% D_bar
  
  # compute F-statistics
  scale_p = (p * (n - 1)) / (n - p)
  F1 = qf(1 - (a), p, n-p)
  F_stat = scale_p * F1
  
  # hypothesis testing
  if (H_T2 > F_stat) {
    cat("Reject null hypothesis\n\n")
    cat("Test-stat:", H_T2, "\n")
    cat("F-quantile:", F_stat)
  } else {
    print("Do not reject null")
  }
  return(data.frame(Hotellings_T2=H_T2, F_qnt=F_stat))
}


## -------------------------------------------------------------------------------------------------------
# compute Bonferroni CI
bonf_CI = function(gr1, gr2, a=.05) {
  # compute parameters
  D = gr1 - gr2
  D_bar = as.matrix(colMeans(D))
  S = cov(D)
  n = dim(gr1)[1]
  p = dim(gr1)[2]
  bon_CI = matrix(rep(0, p * 3), nrow=p)
  for (i in 1:p) {
    bon_CI[i, 1] <- D_bar[i]
    bon_CI[i, 2] <- D_bar[i] - qt(1 - a / (2 * p), n - 1) * sqrt(S[i, i]/n)
    bon_CI[i, 3] <- D_bar[i] + qt(1 - a / (2 * p), n - 1) * sqrt(S[i, i]/n)
  }
  return(data.frame(estimate=bon_CI[, 1], L=bon_CI[, 2], U=bon_CI[, 3]))
}



## -------------------------------------------------------------------------------------------------------
# compare injured to non-injured by paired comparisons
# group 1 (operated)
gr1_I = D_t1_gr1[ , c(2:5)] # injured
gr1_NI = D_t1_gr1[ , c(6:9)] # non-injured
D = gr1_I - gr1_NI


## -------------------------------------------------------------------------------------------------------
m_vec = mahalanobis(D, colMeans(D), cov(D))
m_p = pchisq(m_vec, df=3, lower.tail=FALSE)
idx_rm = which(m_p < 0.05)
m_col = rep("1", 31)
m_col[idx_rm] = "2"
plot(
  1:length(m_p), m_p, ylab="p-value",
  xlab="Observation", col=m_col, main="Outliers in Treatment Group 1"
)
abline(h=0.05, lty=2)

## -------------------------------------------------------------------------------------------------------
c_n = c()
I_n = colnames(gr1_I)
NI_n = colnames(gr1_NI)
for (i in 1:4) {
  c_n[i] = paste(I_n[i], "-", NI_n[i])
}
boxplot(D[-idx_rm, ], main="Treatment Group 1", names=c_n, cex.axis=.8)


## -------------------------------------------------------------------------------------------------------
# check normality group 1
png(file="plot/qq_trt1.png")
par(mfrow=c(2, 2))
for (i in 1:4) {
  qqnorm(
    D[-idx_rm, i], pch=1, frame=FALSE, xlab=if (i > 2) "Theoretical Quantiles" else "",
    ylab=if (i %in% c(1, 3)) "Sample Quantiles" else "",
    main=paste(colnames(gr1_I)[i], "-", colnames(gr1_NI)[i])
  )
  qqline(D[-idx_rm, i], col="3", lwd=2)
  grid()
}


## -------------------------------------------------------------------------------------------------------
# compute test statistic within trt group 1
stat_g1 = paired_comp(gr1_I[-idx_rm,], gr1_NI[-idx_rm, ])


## -------------------------------------------------------------------------------------------------------
# compute CI within trt group 1
CI_g1 = bonf_CI(gr1_I[-idx_rm, ], gr1_NI[-idx_rm, ])
CI_g1

# sign. difference between the test 1 and 2:
# reason to believe that NI > I in these instances

## -------------------------------------------------------------------------------------------------------
print(latex(CI_g1 %>% mutate_if(is.numeric, round, digits=4), file=""))


## -------------------------------------------------------------------------------------------------------
# compare injured to non-injured by paired comparisons
# group 2 (treated but not operated)
gr2_I = D_t1_gr2[ , c(2:5)] # injured
gr2_NI = D_t1_gr2[ , c(6:9)] # non-injured
D = gr2_I - gr2_NI


## -------------------------------------------------------------------------------------------------------
m_vec = mahalanobis(D, colMeans(D), cov(D))
m_p = pchisq(m_vec, df=3, lower.tail=FALSE)
idx_rm = which(m_p < 0.05)
m_col = rep("1", 35)
m_col[idx_rm] = "2"
plot(
  1:length(m_p), m_p, ylab="p-value",
  xlab="Observation", col=m_col, main="Outliers in Treatment Group 2"
)
abline(h=0.05, lty=2)

## -------------------------------------------------------------------------------------------------------
c_n = c()
I_n = colnames(gr1_I)
NI_n = colnames(gr1_NI)
for (i in 1:4) {
  c_n[i] = paste(I_n[i], "-", NI_n[i])
}
boxplot(D[-idx_rm, ], main="Treatment Group 2", names=c_n, cex.axis=.8)


## -------------------------------------------------------------------------------------------------------
# check normality group 2
png(file="plot/qq_trt2.png")
par(mfrow=c(2, 2))
for (i in 1:4) {
  qqnorm(
    D[-idx_rm, i], pch=1, frame=FALSE, xlab=if (i > 2) "Theoretical Quantiles" else "",
    ylab=if (i %in% c(1, 3)) "Sample Quantiles" else "",
    main=paste(colnames(gr1_I)[i], "-", colnames(gr1_NI)[i])
  )
  qqline(D[-idx_rm, i], col="3", lwd=2)
  grid()
}


## -------------------------------------------------------------------------------------------------------
# compute test statistic within trt group 2
stat_g2 = paired_comp(gr2_I[-idx_rm, ], gr2_NI[-idx_rm, ])


## -------------------------------------------------------------------------------------------------------
# compute CI within trt group 2
CI_g2 = bonf_CI(gr2_I[-idx_rm, ], gr2_NI[-idx_rm, ])
CI_g2

# sign. difference between the test 2:
# reason to believe that NI > I in these instances

## -------------------------------------------------------------------------------------------------------
print(latex(CI_g2 %>% mutate_if(is.numeric, round, digits=4), file=""))