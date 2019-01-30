###############################################################################
# Script for pre-processing analysis

library(tidyverse)
library("corrplot")


# Load data
load("rda/galaxy-wd.rda")

# Select only some column of interest for this study
column.select <- c("Type", "D", "D_u", "M_V_T", "M_V_T_u", "N_GC", "N_GC_u",
                   "lg_M_d", "lg_M_d_u", "lg_M_G", "lg_M_G_u")

galaxy <- galaxy %>% select(column.select)

# Relevel galaxies' type
galaxy$Type[grepl("E", galaxy$Type, fixed=FALSE)] <- "E"
galaxy$Type[grepl("S0", galaxy$Type, fixed=FALSE)] <- "S0"
galaxy$Type[grepl("SB0", galaxy$Type, fixed=FALSE)] <- "S0"
galaxy$Type<- droplevels(galaxy$Type)

galaxy$Type[grepl("I", galaxy$Type, fixed=FALSE)] <- "I"
galaxy$Type<- droplevels(galaxy$Type)
not_spiral <- c("E", "I", "S0")
galaxy$Type[ifelse((galaxy$Type %in% not_spiral), FALSE, TRUE)] <- "S"

galaxy$Type<- droplevels(galaxy$Type)

nlevels(galaxy$Type)
table(galaxy$Type)
head(galaxy)

# Compute the intervals
# Some bounds are lower than zero
lower <- galaxy$N_GC-galaxy$N_GC_u
galaxy <- galaxy %>%
  mutate(N_GC.obs.low = ifelse((lower<0.1), 0.1, lower))
galaxy <- galaxy %>%
  mutate(N_GC.obs.up = N_GC + N_GC_u)

# Analyse Data

summary(galaxy)

# Log-Linear plot of globular cluster population

galaxy %>% 
  filter(N_GC > 0) %>%
  ggplot(aes(M_V_T, N_GC)) +
  geom_point(aes(color=Type), alpha=0.6) +
  geom_errorbar(aes(ymin=N_GC-N_GC_u, ymax=N_GC+N_GC_u, 
                    color=Type), width=.2, position=position_dodge(0.05)) +
  geom_errorbarh(aes(xmin=M_V_T-M_V_T_u, xmax=M_V_T+M_V_T_u,
                     color=Type), height=.2) +
  geom_smooth(method = "loess", span = 0.1) +
  scale_x_reverse() +
  scale_y_continuous(trans='log10')+
  #labs(title="Number of Global Clusters", subtitle = "with error bars") +
  labs(y="Globular Clusters") +
  labs(x="Absolute Magnitude") +
  theme(legend.position = "bottom")

ggsave("figs/ClusterLogPlot.pdf", plot = last_plot())


# Global density ditribution of clusters
galaxy %>%
  ggplot(aes(N_GC, fill="g"), color = "grey", show.legend=FALSE) +
  geom_density(alpha = 0.4, show.legend=FALSE) +
  scale_x_continuous(trans='log10') +
  #labs(title="Density ditribution of clusters") +
  labs(y="Density") +
  labs(x="Globular Clusters")

ggsave("figs/PriorDensityMono.pdf", plot = last_plot())

# Density ditribution of clusters given galaxies' types
galaxy %>%
  ggplot(aes(N_GC, fill=Type)) +
  geom_density(alpha = 0.2) +
  scale_x_continuous(trans='log10') +
  #labs(title="Density ditribution of clusters") +
  labs(y="Density") +
  labs(x="Globular Clusters")

ggsave("figs/PriorDensity.pdf", plot = last_plot())

# Global density ditribution of absolute magnitude
# galaxy %>% 
#   ggplot(aes(M_V.T, fill=Type)) + 
#   geom_density(alpha = 0.2) +
#   scale_x_continuous()

# Relevel galaxies' type
galaxy$Type[galaxy$Type=="S0"] <- "S"
galaxy$Type<- droplevels(galaxy$Type)

# Correlation plot
column.select <- c("N_GC", "M_V_T", "D", "lg_M_G")

galaxy %>%
  select(column.select) %>%
  na.omit() %>%
  cor() %>%
  corrplot(type="upper", method="ellipse", tl.pos="d") %>%
  corrplot(type="lower", method="number", col="black", 
           add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")

ggsave("figs/CorrPlot1.pdf", plot = last_plot())

# 2nd correlation plot with fewer columns
column.select <- c("N_GC", "M_V_T", "D")

galaxy %>%
  select(column.select) %>%
  na.omit() %>%
  cor() %>%
  corrplot(type="upper", method="ellipse", tl.pos="d") %>%
  corrplot(type="lower", method="number", col="black", 
           add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")
ggsave("figs/CorrPlot2.pdf", plot = last_plot())

save(galaxy, file = "rda/galaxy.rda")
