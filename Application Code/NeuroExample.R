rm(list = ls())
library(tidyverse)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(gridExtra)
library(grid)
library(ggridges)
library(ggplot2)
library(R.matlab)

datraw <- dplyr::as_tibble(readMat("~/Desktop/Research/Rodu/NPGC/Application/dat.mat")[[1]])
colnames(datraw) <- c("Time", "Stimulus", paste0("Trial ", seq(1:20)))

dat <- datraw %>% filter(Time >= 11.60 & Time < 12.30) %>%
  pivot_longer(cols = 3:22, names_to = "Trial", values_to = "Response")
dat$Trial <- factor(dat$Trial, levels = paste0("Trial ", seq(1:20)))

p1 <- ggplot(dat, aes(x = Time, y = Trial, height = Response + abs(min(Response)))) +
  geom_density_ridges(fill = "white", alpha = 0, stat = "identity") +
  scale_fill_manual(values = c(NA)) +
  scale_x_continuous("Experiment Time (s)", expand = c(0, 0), breaks = seq(11.65, 12.3, 0.1)) +
  scale_y_discrete(paste0("Experimental Trials\nNo. 1 \u2192 20"), expand = c(0.005, 0.005)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 20))

p2 <- ggplot(dat %>% filter(Trial == "Trial 1"), aes(x = Time, y = Stimulus)) +
  geom_line() +
  scale_x_continuous("Experiment Time (s)", expand = c(0, 0), breaks = seq(11.65, 12.3, 0.1)) +
  scale_y_continuous("Signal of\nStimulus", expand = c(0.005, 0.005)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_text(size = 20))

quartz(file = "/Users/noahgade/Desktop/Research/Rodu/NPGC/Figures/app1.pdf", type = "pdf", width = 8, height = 8)
grid.arrange(p2, p1, heights = c(1, 5))
dev.off()

# Example investigating causal delay between signal of stimulus and response in the auditory cortex
# Use 20ms lag for response variable and all predictors
# Select experiment time between 11.60 and 12.30 seconds
sf <- 4000 # sampling frequency (Hz)
# Delays for predictor variables (ms)
delays <- matrix(c(seq(0, 80, 10), seq(20, 100, 10)), ncol = 2)


plot_data <- dplyr::tibble(Delay = seq(1, 9), 
                           Xmin = seq(0, 80, 10),
                           Xmax = seq(20, 100, 10),
                           Ymin = seq(0.8, 0, -0.1),
                           Ymax = seq(0.9, 0.1, -0.1))

plot_data$Value <- c(0.4, 0.18, 0.06, 0.01, 0.01, 0.08, 0.40, 0.74, 0.24)

app2 <- ggplot(plot_data) +
  geom_rect(aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymax, fill = Value), col = "black", size = 0.1) +
  scale_x_continuous("Stimulus-Response Lag (ms)", expand = c(0.02, 0.02), breaks = seq(0, 100, 10)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn("Permutation\nQuantile", 
                      breaks = seq(0, 0.2, 0.04), limits = c(0, 0.2), 
                      colors = c("#CC0000", "#FF9933", rep("#FFFF99", 1)),
                      oob = scales::squish,
                      labels = c("0.00", "0.04", "0.08", "0.12", "0.16", "0.20+")) +
  geom_segment(aes(x = 0, xend = 100, y = 0, yend = 0), size = 0.5) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.key.width = unit(3, "cm"))

quartz(file = "/Users/noahgade/Desktop/Research/Rodu/NPGC/Figures/app2.pdf", type = "pdf", width = 9, height = 4)
app2
dev.off()
