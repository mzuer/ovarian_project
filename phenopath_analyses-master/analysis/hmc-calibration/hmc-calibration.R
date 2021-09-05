
source("simulate_from_phenopath.R")
library(ggplot2)

# Simulate data -----------------------------------------------------------

set.seed(1991)
N <- 100

x <- rep(c(1, -1), each = N / 2)
pst <- rnorm(N)
pst <- scale_vec(pst)

G <- 40

## Let's construct the parameter matrix

de_pars <- t(replicate(G / 4, sample_de()))
pst_pars <- t(replicate(G / 4, sample_pst()))
pst_beta_pars <- t(replicate(G / 4, sample_pst_beta()))
de_pst_beta_pars <- t(replicate(G / 4, sample_de_pst_beta()))

regimes <- c("de", "pst", "pst_beta", "de_pst_beta")

pars <- as_data_frame(rbind(de_pars, pst_pars, pst_beta_pars, de_pst_beta_pars))
names(pars) <- c("alpha", "c", "beta")
pars <- mutate(pars, regime = rep(regimes, each = G / 4))

pars$alpha <- scale_vec(pars$alpha)
pars$c <- scale_vec(pars$c)
pars$beta <- scale_vec(pars$beta)

gex <- apply(pars, 1, function(p) {
  p <- as.numeric(p)
  simulate_one_gene(N, pst, x, p[1], p[2], p[3])
})


# Visualisation -----------------------------------------------------------

library(tidyr)
library(ggplot2)

gex_df <- as_data_frame(gex) %>% 
  mutate(pst, x) %>% 
  gather(feature, expression, -pst, -x)

pca_df <- as_data_frame(prcomp(gex)$x[,1:2]) %>% 
  mutate(pst, x) 

ggplot(pca_df, aes(x = PC1, y = PC2, color = pst)) +
  geom_point() 

ggplot(pca_df, aes(x = PC1, y = PC2, color = factor(x))) +
  geom_point() 


# Inference with CLVM ----------------------------------------------------

library(devtools)
yex <- stat_function(fun = function(x) x, color = 'red')


x <- (x - mean(x)) / sd(x)
pc1 <- pca_df$PC1
pc1 <- (pc1 - mean(pc1)) / sd(pc1)
if(cor(pc1, pst) < 0) pc1 = -pc1

load_all("../../../clvm/")

true_chi <- cbind(rep(10, G), rep(1, G))
true_tau <- rep(1e6, G)

fit <- clvm(scale(gex, scale = F), matrix( 1* x ), pst_init = pst,
            # true_c = pars$c,
            # true_chi = true_chi,
            # true_tau = true_tau,
            # true_beta = matrix(pars$beta, nrow = 1),
            # true_alpha = matrix(pars$alpha, nrow = 1),
            model_mu = F, elbo_tol = 1e-8, maxiter = 10,
            thin = 1, tau_q = 1, tau_c = 1, a = 2, b = 2,
            scale_y = FALSE)#, update_pst = FALSE)
pcavi <- fit

bsf <- fit$beta_sum
bs <- t(fit$m_beta) %*% t(matrix(x))
bs2 <- calculate_greek_sum(fit$m_beta, matrix(x))

bsf[1:3,1:3]
bs[1:3,1:3]
bs2[1:3, 1:3]

asf <- fit$alpha_sum
as <- t(fit$m_alpha) %*% t(matrix(x))
as2 <- calculate_greek_sum(fit$m_alpha, matrix(x))

asf[1:3,1:3]
as[1:3,1:3]
as2[1:3, 1:3]


its <- 1:length(pcavi$elbos)

df <- data_frame(pcavi$elys[its], 
                 pcavi$elps[its], 
                 pcavi$elqs[its], 
                 pcavi$elbos[its], iter = its)#iter = 1:length(pcavi$elys))
names(df) <- c("E[p(y|T)]", "E[p(T)]", "H[q(T)]", "ELBO", "Iter")
gather(df, quantity, value, -Iter) %>%
  ggplot(aes(x = Iter, y = value, color = quantity)) +
  geom_line() +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(subtitle = "Scaled ELBO quantities") + facet_wrap(~ quantity, scales = "free_y")

df <- bind_cols(pcavi[c('ts', 'cs', 'taus', 'alphas', 'betas', 'chis')])
# df <- tbl_df(apply(df, 2, scale_vec))
df$Iter <- 1:nrow(df)
df <- df[its,]

gather(df, quantity, value, -Iter) %>% 
  ggplot(aes(x = Iter, y = value, color = quantity)) +
  geom_line() + # geom_vline(xintercept = 1) +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(subtitle = "Sum of squares quantities") + facet_wrap(~ quantity, scales = "free_y")




plot_elbo(fit)
qplot(pars$alpha, fit$m_alpha[1,])
qplot(pars$alpha, pars$beta - fit$m_beta[1,])


plot_elbo(fit)
elbo <- fit$elbos[-(1:10)]
qplot(seq_along(elbo), elbo, geom = c("point", "line")) + xlab("Iter") + ylab("ELBO")

qplot(pst, fit$m_t, color = factor(x)) + yex
qplot(pars$alpha, fit$m_alpha[1,]) + yex
qplot(pars$c, fit$m_c) + yex
qplot(pars$beta, fit$m_beta[1,]) #+ yex
qplot(x, fit$m_t - pst)
plot(pars$beta - fit$m_beta[1,])

plot(pars$c, fit$m_c - pars$c)
plot(fit$m_c - pars$c)

#qplot(fit$a_tau, fit$b_tau)
qqnorm(fit$m_t)

tau_map <- 1 / (fit$a_tau / fit$b_tau)

g <- which.max(fit$m_beta[1,])

load_all("../../../clvm/")
gs <- calculate_greek_sum(fit$m_alpha, matrix(x))
chi_mat <- matrix(1, nrow = 1, ncol = G)
cavi_update_beta(gs, 0, g - 1, scale(gex, scale = FALSE), matrix(x), fit$m_t, fit$s_t,
                 fit$m_c, fit$m_alpha, fit$m_beta, fit$a_tau,
                 fit$b_tau, chi_mat, chi_mat, rep(0, G))


# Inference with HMC ------------------------------------------------------

library(rstan)

model <- stan_model("~/oxford/cancer/phenopath_analyses/analysis/hmc-calibration/clvm.stan")

## Initialisations
pc1 <- pca_df$PC1
pc1 <- (pc1 - mean(pc1)) / sd(pc1)

c_inits <- apply(gex, 2, function(y) {
  f <- lm(y ~ pc1)
  coef(f)[2]
})
c_inits <- (c_inits - mean(c_inits)) / sd(c_inits)
#t(scale(gex, scale = F))
data <- list(y = t(scale(gex, scale = F)), x = as.vector(x),
             N = nrow(gex), G = ncol(gex),
             z = pst, alpha = pars$alpha, c = pars$c,
             beta = pars$beta, tau = true_tau)

init <- list(list(z = pst,
                 c = c_inits,
                 alpha = rep(0, data$G),
                 beta = rep(0, data$G)))
#sfit <- sampling(model,
#                 data = data, init = init,
#                 chains = 1, iter = 2000)

sfit <- vb(model, data = data, init = init[[1]], grad_samples = 3,
           tol_rel_obj = 0.0005)

poi <- c("beta", "chi", "z", "c", "alpha")#, "tau") # parameters of interest
map_values <- lapply(poi, 
                     function(param) {
                       colMeans(extract(sfit, param)[[param]])
                     })
names(map_values) <- poi

qplot(pars$alpha, pars$alpha - map_values$alpha) 
qplot(pars$beta, map_values$beta) + yex
qplot(pars$beta, pars$beta - map_values$beta) 

qplot(pars$c, map_values$c) + yex
qplot(pars$c, map_values$c - pars$c)

qplot(pst, map_values$z) + yex
qplot(x, map_values$z - pst)#  + yex

qplot(map_values$z, fit$m_t) + yex

qplot(fit$a_tau / fit$b_tau, map_values$tau)



# For scesets -------------------------------------------------------------



library(rstan); set.seed(123)

model <- stan_model("~/oxford/cancer/phenopath_analyses/analysis/hmc-calibration/clvm.stan")

## Initialisations
gex <- scale(t(exprs(sce)))
x <- sce$x
pc1 <- prcomp(gex)$x[,1]
pc1 <- (pc1 - mean(pc1)) / sd(pc1)

c_inits <- apply(gex, 2, function(y) {
  f <- lm(y ~ pc1)
  coef(f)[2]
})
#t(scale(gex, scale = F))
data <- list(y = t(gex), x = as.vector(x),
             N = nrow(gex), G = ncol(gex))
init <- list(list(z = pc1,
                  c = c_inits,
                  alpha = rep(0, data$G),
                  beta = rep(0, data$G)))

sfit <- sampling(model,
                 data = data, init = init,
                 chains = 1, iter = 2000)

sfit <- vb(model, data = data, init = init[[1]], grad_samples = 1)

poi <- c("alpha", "beta", "c", "z", "chi", "tau") # parameters of interest
map_values <- lapply(poi, 
                     function(param) {
                       colMeans(extract(sfit, param)[[param]])
                     })
names(map_values) <- poi

qplot(pc1, map_values$z, color = x)

fit <- clvm(gex, matrix(x), pst_init = pc1, 
            model_mu = FALSE, elbo_tol = 1e-8)

elbo <- fit$elbos[-(1)]
qplot(seq_along(elbo) + 1, elbo, geom = c("point", "line")) + xlab("Iter") + ylab("ELBO")

load_all("../../../clvm/")
pc1u <- prcomp(gex)$x[,1]
fit2 <- clvm(gex, matrix(x), pst_init = pc1u, tau_c = 10, 
             #true_beta = matrix(map_values$beta, nrow =1),
             #true_c = map_values$c,
            model_mu = FALSE, elbo_tol = 1e-3)

qplot(map_values$z, fit$m_t) +yex +# stat_function(fun = function(x) x, color = 'red') +
  #stat_function(fun = function(x) 2 * x, color = 'green') +
  xlab("Stan") + ylab("Cavi")

qplot(map_values$c, fit$m_c) + yex
qplot(map_values$alpha, fit$m_alpha[1,]) + yex
qplot(map_values$beta, fit$m_beta[1,]) + yex

gs <- clvm:::calculate_greek_sum(fit$m_beta, matrix(x))


# Beta by hand ------------------------------------------------------------
yscaled <- scale(gex, scale = FALSE)

sample_t_sq <- function(i) {
  rnorm(500, fit$m_t[i], sqrt(fit$s_t[i]))^2
}

tsq_mc <- sapply(1:N, function(i) mean(sample_t_sq(i) * x[i]^2))

# plot(tsq_mc, )

s_beta <- fit$a_tau / fit$b_tau * sum(fit$m_t^2 + fit$s_t) * x^2 + fit$chi_exp[1,]

m_beta <- sapply(1:G, function(g) {
  r <- 0
  for(i in 1:N) {
    r <- r + fit$m_t[i] * x[i] * yscaled[i,g] -
      x[i] * fit$m_c[g] * (fit$m_t[i]^2 + fit$s_t[i]) -
      fit$m_alpha[1,g] * x[i]^2 * fit$m_t[i]
  }
  r <- r * fit$a_tau[g] / fit$b_tau[g]
  r
})

m_beta <- m_beta / s_beta

alpha_sum <- clvm:::calculate_greek_sum(fit$m_alpha, matrix(x))
devtools::load_all("../../../clvm/")

z <- clvm:::cavi_update_beta(alpha_sum, 0, 0, yscaled, matrix(x),
                        fit$m_t, fit$s_t, fit$m_c, fit$m_alpha,
                        fit$m_beta, fit$a_tau, fit$b_tau,
                        t(true_chi[,1,drop=F]), t(true_chi[,2,drop=F]),
                        fit$m_mu)
z[2]
mmb <- z[1] / z[2]

# Beta expectation --------------------------------------------------------

be <- be2 <- matrix(0, nrow = N, ncol = G)

yscaled <- scale(gex, scale = FALSE)
for(i in 1:N) {
  for(g in 1:G) {
    tmp <- fit$m_c[g]^2 + fit$s_c[g]
    tmp <- tmp + 2 * fit$m_c[g] * fit$m_beta[1,g] * x[i]
    tmp <- tmp + (fit$m_beta[1,g]^2 + fit$s_beta[1,g]) * x[i]^2
    be[i,g] <- tmp
  }
}


# MC it
sample_beta <- function(g) rnorm(1000, fit$m_beta[1,g], sqrt(fit$s_beta[1,g]))
sample_c <- function(g) rnorm(1000, fit$m_c[g], sqrt(fit$s_c[g]))

for(i in 1:N) {
  for(g in 1:G) {
    b <- sample_beta(g) * x[i] + sample_c(g) 
    be2[i,g] <- mean(b^2)
  }
}

qplot(factor(x), rowMeans(be2), geom = 'box') +
  ylab(expression(Sigma[g]~E~"["~c[g]~+~beta[g]~"]"))
  
qplot(as.vector(be), as.vector(be2)) + yex

## MC for numerator
sample_beta <- function(g) rnorm(500, fit$m_beta[1,g], sqrt(fit$s_beta[1,g]))
sample_c <- function(g) rnorm(500, fit$m_c[g], sqrt(fit$s_c[g]))
sample_alpha <- function(g) rnorm(500, fit$m_alpha[1,g], sqrt(fit$s_alpha[1,]))

num <- matrix(0, nrow = N, ncol = G)
for(i in 1:N) {
  for(g in 1:G) {
    n <- (sample_c(g) + sample_beta(g) * x[i] ) * (yscaled[i,g] - sample_alpha(g))
    num[i,g] <- mean(n)
  }
}

qplot(factor(x), rowMeans(num), geom = 'boxplot')



m_t <- rep(0, N)
s_t <- rep(1, N)
yscaled <- scale(gex, scale = FALSE)
for(i in 1:N) {
  for(g in 1:G) {
    tmp <- fit$m_c[g]^2 + fit$s_c[g]
    tmp <- tmp + 2 * fit$m_c[g] * fit$m_beta[1,g] * x[i]
    tmp <- tmp + (fit$m_beta[1,g]^2 + fit$s_beta[1,g]) * x[i]^2
    s_t[i] <- s_t[i] + fit$a_tau[g] / fit$b_tau[g] * tmp
    
    m <- (yscaled[i,g] - fit$m_alpha[1,g] * x[i]) * (fit$m_c[g] + fit$m_beta[1,g] * x[i])
    m_t[i] <- m_t[i] + fit$a_tau[g] / fit$b_tau[g] * m
  }
}

m_t <- m_t / s_t
s_t <- 1 / s_t

pup <- cavi_update_pst(scale(gex, scale = F), matrix(x), fit$m_c, fit$m_mu, fit$s_c,
                       fit$m_alpha, fit$m_beta, fit$s_beta,
                       fit$a_tau, fit$b_tau, rep(0, N), 1)

plot(x, m_t - pst)
plot(x, fit$m_t - pst)
plot(x, pup[,1] - pst)

Ebx <- clvm:::calculate_greek_sum(fit$m_beta, matrix(x))
Ebx2 <- clvm:::greek_square_exp(fit$m_beta, fit$s_beta, matrix(x))

l <- rep(0, N)
k <- rep(0, N)
for(i in 1:N) {
  for(g in 1:G) {
    l[i] <- l[i] + (gex[i,g] - fit$m_alpha[1,g] * x[i])
    k[i] <- k[i] + fit$m_c[g] + fit$m_beta[1,g] * x[i]  
  }
}




# ELBO --------------------------------------------------------------------
library(tidyverse)
scale_vec <- function(x)  (x - mean(x)) /sd(x)
df <- data_frame(scale_vec(fit$elys), 
                 scale_vec(fit$elps), 
                 scale_vec(fit$elqs), iter = 1:length(fit$elys))
names(df) <- c("E[p(y|T)]", "E[p(T)]", "H[q(T)]", "Iter")
df <- df[1:20,]
gather(df, quantity, value, -Iter) %>% 
  ggplot(aes(x = Iter, y = value, color = quantity)) +
  geom_line() +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(subtitle = "Scaled ELBO quantities")


# Check fg calculation ----------------------------------------------------

sample_beta <- function(g) rnorm(1, fit$m_beta[1,g], sqrt(fit$s_beta[1,g]))
sample_c <- function(g) rnorm(1, fit$m_c[g], sqrt(fit$s_c[g]))
sample_alpha <- function(g) rnorm(1, fit$m_alpha[1,g], sqrt(fit$s_alpha[1,]))
sample_t <- function() rnorm(length(fit$m_t), fit$m_t, sqrt(fit$s_t))
sample_tau <- function(g) rgamma(1, fit$a_tau[g], fit$b_tau[g])


alpha_sum <- clvm:::calculate_greek_sum(fit$m_alpha, matrix(x))
alpha_square_sum <- clvm:::greek_square_exp(fit$m_alpha, fit$s_alpha, matrix(x))

beta_sum <- clvm:::calculate_greek_sum(fit$m_beta, matrix(x))
beta_square_sum <- clvm:::greek_square_exp(fit$m_beta, fit$s_beta, matrix(x))

g <- 56
fg_analytical <- calculate_fg(g-1, yscaled, fit$m_t, fit$s_t, fit$m_c, fit$s_c, fit$m_alpha, fit$s_alpha, fit$m_beta, 
                              fit$s_beta, rep(0,G), rep(0,G), 
                              alpha_sum, beta_sum, alpha_square_sum, beta_square_sum)

f <- sapply(1:G, function(g) {
  fit$a_tau[g] / (2 * fit$b_tau[g]) * calculate_fg(g-1, yscaled, fit$m_t, fit$s_t, fit$m_c, fit$s_c, fit$m_alpha, fit$s_alpha, fit$m_beta, 
               fit$s_beta, rep(0,G), rep(0,G), 
               alpha_sum, beta_sum, alpha_square_sum, beta_square_sum)
})

f_g <- function(y, alpha, beta, x, c, t) {
  sum((y - alpha * x - c * t - beta * x * t)^2)
}

fg_devs <- fg_analytical - replicate(100000,
          f_g(yscaled[,g], sample_alpha(g), sample_beta(g), 
              x, sample_c(g), sample_t()))
qplot(fg_devs, geom = 'density')

calculate_E_log_Y_given_theta(yscaled, matrix(x), fit$m_t, fit$s_t, fit$m_c,
                              fit$s_c, fit$m_alpha, fit$s_alpha, fit$m_beta, 
                              fit$s_beta, fit$a_tau, fit$b_tau,
                              rep(0,G), rep(0,G))

mu <- function(alpha, beta, x, c, t) {
  alpha * x + c * t + beta * x * t
}

sample_t2 <- function(i) rnorm(1, fit$m_t[i], sqrt(fit$s_t[i]))

log_p_y <- function() {
  lpy <- 0
  for(i in 1:N) {
    for(g in 1:G) {
      m <- mu(sample_alpha(g), sample_beta(g),
              x[i], sample_c(g), sample_t2(i))
      prec <- sample_tau(g)
      lpy <- lpy + dnorm(yscaled[i,g], m, 1 / sqrt(prec), log = TRUE)
    }
  }
  return(lpy)
}

lpy_mc <- mean(replicate(250, log_p_y()))
lpy_mc

elt_mc <- replicate(1000, sample_tau(1))
mean(log(elt_mc))

digamma(fit$a_tau[1]) - log(fit$b_tau[1])

v0 <- rep(0, N)
m0 <- matrix(0, nrow = 1, ncol = G)
m1 <- matrix(0, nrow = N, ncol = G)

ys <- matrix(1)

# g, y, m_t, s_t, m_c, s_c, m_alpha, s_alpha, 
# m_beta, s_beta, m_mu, s_mu, alpha_sum, beta_sum, 
# alpha_square_sum, beta_square_sum)
calculate_fg(0, ys, 1, 1, 1, 1,
             m0, m0, m0, m0, v0, v0,
            m1, m1, m1, m1)

zzz <- replicate(100000,
          f_g(1, 0, 0, 
              0, rnorm(1, 1, 1), rnorm(1, 1, 1)))


#### Take this down a level
m_alpha <- 0
s_alpha <- 1
m_beta <- 0
s_beta <- 1
m_t <- 0
s_t <- 1
y <- 1
m_c <- 1
s_c <- 1
xx <- 1
m_mu <- 0
s_mu <- 0
a <- 2; b <- 2
alpha_sum <- m_alpha * xx; beta_sum <- m_beta * xx
alpha_square_sum <- xx^2 * (m_alpha^2 + s_alpha)
beta_square_sum <- xx^2 * (m_beta^2 + s_beta)
mu2 <- function() {
  t <- rnorm(1, m_t, sqrt(s_t))
  rnorm(1, m_alpha, sqrt(s_alpha)) * xx + 
    rnorm(1, m_c, sqrt(s_c)) * t + 
    rnorm(1, m_beta, sqrt(s_beta)) * xx * t
}
ff <- function() (y - mu2())^2


# g, y, m_t, s_t, m_c, s_c, m_alpha, s_alpha, 
# m_beta, s_beta, m_mu, s_mu, alpha_sum, beta_sum, 
# alpha_square_sum, beta_square_sum)
calculate_fg(0, matrix(y), m_t, s_t, m_c, s_c, 0, 0, 
             matrix(alpha_sum), matrix(beta_sum), 
             matrix(alpha_square_sum), matrix(beta_square_sum))

fs <- replicate(1e6, ff())
mean(fs)

log_p_y <- function() {
  lpy <- 0
  m <- mu2()
  prec <- rgamma(1, 2, 1)
  dnorm(y, m, 1 / sqrt(prec), log = TRUE)
}
calculate_E_log_Y_given_theta(matrix(y), matrix(xx), m_t, s_t, m_c,
                              s_c, matrix(m_alpha), matrix(s_alpha), matrix(m_beta),
                              matrix(s_beta), 2, 1,
                              0, 0)
lpys <- replicate(5e5, log_p_y())
mean(lpys)


# Beta sum ----------------------------------------------------------------

bs <- t(fit$m_beta) %*% t(matrix(x))
bs2 <- calculate_greek_sum(fit$m_beta, matrix(x))
