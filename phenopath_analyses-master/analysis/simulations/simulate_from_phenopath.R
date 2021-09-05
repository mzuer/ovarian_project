
library(dplyr)
library(polyester)
library(Biostrings)
library(readr)

set.seed(123L)

mutate <- dplyr::mutate

simulate_one_gene <- function(N, pst, x, alpha = 0, c = 0, beta = 0, tau = 1e6) {
  mu <- alpha * x + (c + beta * x) * pst
  rnorm(N, mu, 1 / sqrt(tau))
}

sample_de <- function() {
  alpha <- sample(c(-1, 1), 1)
  c <- beta <- 0
  return(c(alpha, c, beta))
}

sample_pst <- function() {
  c <- 1 * sample(c(-1, 1), 1)
  alpha <- beta <- 0
  return(c(alpha, c, beta))
}

sample_pst_beta <- function() {
  c <- 1 * sample(c(-1, 1), 1)  
  beta <- sample(c(-1, 1), 1)  
  alpha <- 0
  return(c(alpha, c, beta))
}

sample_de_pst_beta <- function() {
  c <- 1 * sample(c(-1, 1), 1)  
  beta <- sample(c(-1, 1), 1)  
  alpha <- sample(c(-1, 1), 1)
  return(c(alpha, c, beta))
}


# Code starts here --------------------------------------------------------


N <- 200

x <- rep(c(1, -1), each = N / 2)
pst <- rnorm(N)

G <- 400

## Let's construct the parameter matrix

de_pars <- t(replicate(G / 4, sample_de()))
pst_pars <- t(replicate(G / 4, sample_pst()))
pst_beta_pars <- t(replicate(G / 4, sample_pst_beta()))
de_pst_beta_pars <- t(replicate(G / 4, sample_de_pst_beta()))

regimes <- c("de", "pst", "pst_beta", "de_pst_beta")

pars <- as_data_frame(rbind(de_pars, pst_pars, pst_beta_pars, de_pst_beta_pars))
names(pars) <- c("alpha", "c", "beta")
pars <- mutate(pars, regime = rep(regimes, each = G / 4))

gex <- apply(pars, 1, function(p) {
  p <- as.numeric(p)
  simulate_one_gene(N, pst, x, p[1], p[2], p[3])
})

# stop("done")

pos_gex <- 2 ^ gex

count_mat <- sapply(seq_len(nrow(pos_gex)), function(i) {
  x <- pos_gex[i, ]
  NB(x, x / 3)
})

write_csv(as_data_frame(count_mat), "../../data/simulations/count_matrix.csv")

# FASTA annotation
fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)

# subset the FASTA file to first 20 transcripts
small_fasta = fasta[1:G]

ref_file <- "data/simulations/ref/chr22_small.fa"

writeXStringSet(small_fasta, ref_file)


simulate_experiment_countmat(ref_file, readmat = count_mat, 
                             outdir = "data/simulations/fasta/")



write_csv(pars, "data/simulations/gene_pars.csv")

pdata <- data_frame(sample = 1:N, pst, x)

write_csv(pdata, "data/simulations/pdata.csv")

print("Done simulating experiment")
