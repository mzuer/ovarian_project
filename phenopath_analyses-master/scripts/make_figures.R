library(tidyverse)
library(viridis)
library(cowplot)

sigmoid <- function(pst, mu0, k, t0) 2 * mu0 / (1 + exp(-k*(pst - t0)))
linear <- function(pst, c, k) c + k*pst

N <- 100

pst <- runif(N)


# Discrete case -----------------------------------------------------------

df_discrete_predicted <- data_frame(
  pst = pst,
  x1 = sigmoid(pst, 2, 10, 0.7),
  x2 = sigmoid(pst, -2, 10, 0.7)
)

df_discrete_noise <- mutate(df_discrete_predicted,
                      x1 = rnorm(length(x1), x1, .6),
                      x2 = rnorm(length(x2), x2, .6))

df_pred <- gather(df_discrete_predicted, covariate, y, -pst) %>% arrange(pst)
df_noise <- gather(df_discrete_noise, covariate, y, -pst)


ggplot(df_pred, aes(x = pst, y = y, color = covariate)) +
  #geom_point(data = df_noise, alpha = 0.4) +
  geom_line(size = 1.5) +
  scale_color_brewer(name = "Discrete\nCovariate", palette = "Set1") +
  xlab("Trajectory") + ylab("Gene expression") +
  theme(legend.text = element_blank(), #legend.direction = "vertical",
        legend.position = "top", axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_text(size = 11),
        legend.title = element_text(size = 11))

discrete_plot <- last_plot()



# Continuous case ---------------------------------------------------------

mu0 <- seq(-2, 2, length.out = 100)

df_conts <- lapply(mu0, function(m0) {
  data_frame(
    pst = pst,
    y = sigmoid(pst, m0, 10, 0.7),
    mu0 = m0
  )
})

dfc <- bind_rows(df_conts)

# dfc_noise <- mutate(dfc, y = rnorm(length(y), y, .6))

dfc %>% arrange(mu0, pst) %>% 
  ggplot(aes(x = pst, y = y, group = mu0, color = mu0)) +
  # geom_point(data = dfc_noise) +
  geom_line(size = 1.5) + 
  scale_color_viridis(name = "Continuous\nCovariate") +
  theme(legend.text = element_blank(), # legend.direction = "horizontal",
        legend.title = element_text(size = 11),
        legend.position = "top", axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_text(size = 11),
        legend.key.size = unit(.2, "cm")) +
  xlab("Trajectory") + ylab("Gene expression")

continuous_plot = last_plot()

plot_grid(discrete_plot, continuous_plot, rel_widths = c(10, 11))

phenotime_diagram <- last_plot()

# save(phenotime_diagram, file = "../figs/phenotime_diagram.Rdata")
ggsave("../figs/method_diagram.png", width = 7, height = 3.5)



# Linear discrete ---------------------------------------------------------

df_discrete <- data_frame(
  pst = pst,
  x1 = linear(pst, 2, -2),
  x2 = linear(pst, 2, 0)
)

gather(df_discrete, covariate, y, -pst) %>% 
  arrange(pst) %>% 
  ggplot(aes(x = pst, y = y, color = covariate)) +
  geom_line(size = 1.5) +
  scale_color_brewer(name = "Discrete\nCovariate (x)", palette = "Set1") +
  xlab("Patient trajectory (z)") + ylab("Dynamic observable (y)") +
  theme(legend.text = element_blank(), legend.direction = "horizontal",
        legend.position = "top", axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_text(size = 11))

discrete_plot <- last_plot()


# Continuous linear -------------------------------------------------------


ks <- seq(-2, 2, length.out = 100)

df_conts <- lapply(ks, function(k) {
  data_frame(
    pst = pst,
    y = linear(pst, 0, k),
    k = k
  )
})

dfc <- bind_rows(df_conts)

dfc %>% arrange(k, pst) %>% 
  ggplot(aes(x = pst, y = y, group = k, color = k)) +
  geom_line(size = 1.5) + scale_color_viridis(name = "Continuous\nCovariate (x)") +
  theme(legend.text = element_blank(), legend.direction = "horizontal",
        legend.title = element_text(hjust = 0.5),
        legend.position = "top", axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_text(size = 11)) +
  xlab("Patient trajectory (z)") + ylab("Dynamic observable (y)")

continuous_plot = last_plot()

plot_grid(discrete_plot, continuous_plot)

phenotime_diagram_linear <- last_plot()

save(phenotime_diagram_linear, file = "../figs/phenotime_diagram_linear.Rdata")
ggsave("../figs/phenotime_diagram_linear.png", width = 9, height = 3)
  
  