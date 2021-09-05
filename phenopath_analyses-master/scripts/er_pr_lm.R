
sce <- readRDS("data/BRCA/sce_brca_clvm.rds")
pcavi <- readRDS("data/BRCA/clvm_results.rds")

fit_model <- function(y, ER, PR, HER2, z) {
  fit <- lm(y ~ ER*z + PR*z + HER2 * z)
  coef(fit)
}

scale_vec <- function(x) (x - mean(x)) / sd(x)

ER <- scale_vec(1 * (sce$ER_status == "positive"))
PR <- scale_vec(1 * (sce$PR_status == "positive"))

HER2 <- sce$IHC_HER2
is_equiv <- HER2 == "Equivocal" | HER2 == "Indeterminate" | is.na(HER2)

h2 <- HER2
#h2[is_equiv] <- 0
#h2[!is_equiv] <- scale_vec(1 * (h2[!is_equiv] == "Positive"))
h2[is.na(h2)] <- "NA"
h2 <- scale_vec(1 * (h2 == "Positive"))

z <- scale_vec(pcavi$m_t)

coefs <- apply(scale(t(exprs(sce))), 2, fit_model, ER, PR, h2, z)

plot(pcavi$m_beta[1,], coefs[5,])

beta_er <- coefs[6,]
beta_pr <- coefs[7,]
beta_her2 <- coefs[8,]

df <- data.frame(beta_er, beta_pr, beta_her2)


qplot(beta_er, beta_pr, alpha = 0.5)



