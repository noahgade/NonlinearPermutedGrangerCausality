rm(list = ls())

library(tidyverse)

# Read in Simulation Data
load("~/Desktop/Research/Rodu/NPGC/Simulation Results/naive250.RData")
load("~/Desktop/Research/Rodu/NPGC/Simulation Results/naive500.RData")
load("~/Desktop/Research/Rodu/NPGC/Simulation Results/naive1000.RData")
load("~/Desktop/Research/Rodu/NPGC/Simulation Results/npgc250.RData")
load("~/Desktop/Research/Rodu/NPGC/Simulation Results/npgc500.RData")
for(i in 1:8) {
  load(paste0("~/Desktop/Research/Rodu/NPGC/Simulation Results/npgc1000_", i, ".RData"))
  assign(paste0("npgc1000_", i), npgc1000)
}
npgc1000 <- npgc1000_1 %>% append(npgc1000_2) %>% append(npgc1000_3) %>% append(npgc1000_4) %>%
  append(npgc1000_5) %>% append(npgc1000_6) %>% append(npgc1000_7) %>% append(npgc1000_8)

GL250 <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/GL250.csv")
GL500 <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/GL500.csv")
GL1000 <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/GL1000.csv")

GSGL250 <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/GSGL250.csv")
GSGL500 <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/GSGL500.csv")
GSGL1000 <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/GSGL1000.csv")

H250 <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/H250.csv")
H500 <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/H500.csv")
H1000 <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/H1000.csv")

extractQres <- function(list_data, len, n_initializations) {
  statistic <- mean(list_data[2,] / list_data[1,])
  m <- (len - 100) / (len - 102)
  s <- sqrt((2 * (len - 100) ^ 2 * (2 * len - 202)) / ((len - 100) * (len - 102) ^ 2 * (len - 104) * n_initializations))
  quant <- 1 - pnorm(statistic, m, s)
  return(quant)
}

extractQnoise <- function(list_data, len, n_initializations) {
  statistic <- mean(list_data[3,] / list_data[1,])
  m <- (len - 100) / (len - 102)
  s <- sqrt((2 * (len - 100) ^ 2 * (2 * len - 202)) / ((len - 100) * (len - 102) ^ 2 * (len - 104) * n_initializations))
  quant <- 1 - pnorm(statistic, m, s)
  return(quant)
}

extractQnpgc <- function(npgc_result) {
  out <- apply(npgc_result, 3, mean)
  o <- sum(out <= out[1]) / length(out)
  return(o)
}

extractCMLP <- function(cmlpresult) {
  out <- apply(cmlpresult, 1, function(x) ifelse(sum(x[4:7]) == 0, 1, 0))
  return(out)
}

# Extracting results by dataset
## T = 250
Results250 <- dplyr::tibble(NPGC = unlist(lapply(npgc250, extractQnpgc)),
                            cmlpGL = extractCMLP(GL250),
                            cmlpGSGL = extractCMLP(GSGL250),
                            cmlpH = extractCMLP(H250),
                            RU = unlist(lapply(naive250, extractQres, len = 250, n_initializations = 1000)),
                            GNS = unlist(lapply(naive250, extractQnoise, len = 250, n_initializations = 1000)),
                            TimePts = 250,
                            Type = rep(c(rep("Lorenz-96", 200), rep("TAR(2)", 200)), 4),
                            Causal = c(rep("NC", 800), rep("GC", 800)),
                            NZ = rep(c(rep(0, 400), rep(2, 400)), 2))

## T = 500
Results500 <- dplyr::tibble(NPGC = unlist(lapply(npgc500, extractQnpgc)),
                            cmlpGL = extractCMLP(GL500),
                            cmlpGSGL = extractCMLP(GSGL500),
                            cmlpH = extractCMLP(H500),
                            RU = unlist(lapply(naive500, extractQres, len = 500, n_initializations = 1000)),
                            GNS = unlist(lapply(naive500, extractQnoise, len = 500, n_initializations = 1000)),
                            TimePts = 500,
                            Type = rep(c(rep("Lorenz-96", 200), rep("TAR(2)", 200)), 4),
                            Causal = c(rep("NC", 800), rep("GC", 800)),
                            NZ = rep(c(rep(0, 400), rep(2, 400)), 2))

## T = 1000
Results1000 <- dplyr::tibble(NPGC = unlist(lapply(npgc1000, extractQnpgc)),
                            cmlpGL = extractCMLP(GL1000),
                            cmlpGSGL = extractCMLP(GSGL1000),
                            cmlpH = extractCMLP(H1000),
                            RU = unlist(lapply(naive1000, extractQres, len = 1000, n_initializations = 1000)),
                            GNS = unlist(lapply(naive1000, extractQnoise, len = 1000, n_initializations = 1000)),
                            TimePts = 1000,
                            Type = rep(c(rep("Lorenz-96", 200), rep("TAR(2)", 200)), 4),
                            Causal = c(rep("NC", 800), rep("GC", 800)),
                            NZ = rep(c(rep(0, 400), rep(2, 400)), 2))

T1error <- dplyr::tibble(Reference = rep(seq(1/400, 399/400, 1/400), 36),
                         Method = rep(rep(c("NPGC", "cmlpGL", "cmlpGSGL", "cmlpH", "RU", "GNS"), each = 399), 6),
                         TimePts = rep(rep(c(250, 500, 1000), each = 2394), 2),
                         Type = c(rep("Lorenz-96", 7182), rep("TAR(2)", 7182)),
                         Observed = NA)

for(type in c("Lorenz-96", "TAR(2)")) {
  for(method in c("NPGC", "cmlpGL", "cmlpGSGL", "cmlpH", "RU", "GNS")) {
    temp250 <- Results250 %>% 
      pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>%
      filter(Type == type & Causal == "NC" & Method == method) %>% select(Q)
    temp500 <- Results500 %>% 
      pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>% 
      filter(Type == type & Causal == "NC" & Method == method) %>% select(Q)
    temp1000 <- Results1000 %>%
    pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>%
     filter(Type == type & Causal == "NC" & Method == method) %>% select(Q)
    
    T1error[(T1error$Type == type & T1error$Method == method & T1error$TimePts == 250), ]$Observed <- sapply(seq(1/400, 399/400, 1/400), function(x) sum(temp250 <= x) / 400)
    T1error[(T1error$Type == type & T1error$Method == method & T1error$TimePts == 500), ]$Observed <- sapply(seq(1/400, 399/400, 1/400), function(x) sum(temp500 <= x) / 400)
    T1error[(T1error$Type == type & T1error$Method == method & T1error$TimePts == 1000), ]$Observed <- sapply(seq(1/400, 399/400, 1/400), function(x) sum(temp1000 <= x) / 400)
  }
}

## Figure 1a
fig1a <- ggplot(T1error %>% filter(Type == "TAR(2)")) +
  geom_line(aes(x = Reference, y = Observed, color = Method, linetype = Method)) +
  scale_color_manual("Method", values = c("NPGC" = "#D62728", "cmlpGL" = "#1F77B4", "cmlpGSGL" = "#1F77B4",
                                          "cmlpH" = "#1F77B4", "RU" = "#FF7F0E", "GNS" = "#2CA02C")) +
  scale_linetype_manual("Method", values = c("NPGC" = "solid", "cmlpGL" = "dotted", "cmlpGSGL" = "dashed",
                                             "cmlpH" = "solid", "RU" = "solid", "GNS" = "solid")) +
  geom_abline(aes(slope = 1, intercept = 0), color = "lightgray") +
  facet_wrap(~ TimePts, nrow = 1,
             labeller = labeller(TimePts = c("250" = "T = 250", "500" = "T = 500", "1000" = "T = 1000"))) +
  scale_x_continuous("Level of Test", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  scale_y_continuous("Observed Type 1 Error", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank())

quartz(file = "/Users/noahgade/Desktop/Research/Rodu/NPGC/Figures/figure1a.pdf", type = "pdf", width = 13, height = 5)
fig1a
dev.off()

fig1b <- ggplot(T1error %>% filter(Type == "Lorenz-96")) +
  geom_line(aes(x = Reference, y = Observed, color = Method, linetype = Method)) +
  scale_color_manual("Method", values = c("NPGC" = "#D62728", "cmlpGL" = "#1F77B4", "cmlpGSGL" = "#1F77B4",
                                          "cmlpH" = "#1F77B4", "RU" = "#FF7F0E", "GNS" = "#2CA02C")) +
  scale_linetype_manual("Method", values = c("NPGC" = "solid", "cmlpGL" = "dotted", "cmlpGSGL" = "dashed",
                                             "cmlpH" = "solid", "RU" = "solid", "GNS" = "solid")) +
  geom_abline(aes(slope = 1, intercept = 0), color = "lightgray") +
  facet_wrap(~ TimePts, nrow = 1,
             labeller = labeller(TimePts = c("250" = "T = 250", "500" = "T = 500", "1000" = "T = 1000"))) +
  scale_x_continuous("Level of Test", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  scale_y_continuous("Observed Type 1 Error", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank())

quartz(file = "/Users/noahgade/Desktop/Research/Rodu/NPGC/Figures/figure1b.pdf", type = "pdf", width = 13, height = 5)
fig1b
dev.off()




T1errorB <- dplyr::tibble(Reference = rep(seq(1/200, 199/200, 1/200), 72),
                         Method = rep(rep(c("NPGC", "cmlpGL", "cmlpGSGL", "cmlpH", "RU", "GNS"), each = 199), 12),
                         TimePts = rep(rep(c(250, 500, 1000), each = 2388), 2),
                         Type = c(rep("Lorenz-96", 7164), rep("TAR(2)", 7164)),
                         NZ = rep(rep(c(0, 2), each = 1194), 6),
                         Observed = NA)

for(type in c("Lorenz-96", "TAR(2)")) {
  for(method in c("NPGC", "cmlpGL", "cmlpGSGL", "cmlpH", "RU", "GNS")) {
    for(nz in c(0,2)) {
      temp250 <- Results250 %>% 
        pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>%
        filter(Type == type & Causal == "NC" & Method == method & NZ == nz) %>% select(Q)
      
      temp500 <- Results500 %>% 
        pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>%
        filter(Type == type & Causal == "NC" & Method == method & NZ == nz) %>% select(Q)
      
      temp1000 <- Results1000 %>% 
        pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>%
        filter(Type == type & Causal == "NC" & Method == method & NZ == nz) %>% select(Q)
      
      T1errorB[(T1errorB$Type == type & T1errorB$Method == method & T1errorB$TimePts == 250 & T1errorB$NZ == nz), ]$Observed <- sapply(seq(1/200, 199/200, 1/200), function(x) sum(temp250 <= x) / 200)
      T1errorB[(T1errorB$Type == type & T1errorB$Method == method & T1errorB$TimePts == 500 & T1errorB$NZ == nz), ]$Observed <- sapply(seq(1/200, 199/200, 1/200), function(x) sum(temp500 <= x) / 200)
      T1errorB[(T1errorB$Type == type & T1errorB$Method == method & T1errorB$TimePts == 1000 & T1errorB$NZ == nz), ]$Observed <- sapply(seq(1/200, 199/200, 1/200), function(x) sum(temp1000 <= x) / 200)
    }
  }
}

figA <- ggplot(T1errorB %>% filter(Type == "TAR(2)" & NZ == 0)) +
  geom_line(aes(x = Reference, y = Observed, color = Method, linetype = Method)) +
  scale_color_manual("Method", values = c("NPGC" = "#D62728", "cmlpGL" = "#1F77B4", "cmlpGSGL" = "#1F77B4",
                                          "cmlpH" = "#1F77B4", "RU" = "#FF7F0E", "GNS" = "#2CA02C")) +
  scale_linetype_manual("Method", values = c("NPGC" = "solid", "cmlpGL" = "dotted", "cmlpGSGL" = "dashed",
                                             "cmlpH" = "solid", "RU" = "solid", "GNS" = "solid")) +
  geom_abline(aes(slope = 1, intercept = 0), color = "lightgray") +
  facet_wrap(~ TimePts, nrow = 1,
             labeller = labeller(TimePts = c("250" = "T = 250", "500" = "T = 500", "1000" = "T = 1000"))) +
  scale_x_continuous("Level of Test", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  scale_y_continuous("Observed Type 1 Error", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank())

quartz(file = "/Users/noahgade/Desktop/Research/Rodu/NPGC/Figures/figureA.pdf", type = "pdf", width = 13, height = 5)
figA
dev.off()

figB <- ggplot(T1errorB %>% filter(Type == "TAR(2)" & NZ == 2)) +
  geom_line(aes(x = Reference, y = Observed, color = Method, linetype = Method)) +
  scale_color_manual("Method", values = c("NPGC" = "#D62728", "cmlpGL" = "#1F77B4", "cmlpGSGL" = "#1F77B4",
                                          "cmlpH" = "#1F77B4", "RU" = "#FF7F0E", "GNS" = "#2CA02C")) +
  scale_linetype_manual("Method", values = c("NPGC" = "solid", "cmlpGL" = "dotted", "cmlpGSGL" = "dashed",
                                             "cmlpH" = "solid", "RU" = "solid", "GNS" = "solid")) +
  geom_abline(aes(slope = 1, intercept = 0), color = "lightgray") +
  facet_wrap(~ TimePts, nrow = 1,
             labeller = labeller(TimePts = c("250" = "T = 250", "500" = "T = 500", "1000" = "T = 1000"))) +
  scale_x_continuous("Level of Test", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  scale_y_continuous("Observed Type 1 Error", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank())

quartz(file = "/Users/noahgade/Desktop/Research/Rodu/NPGC/Figures/figureB.pdf", type = "pdf", width = 13, height = 5)
figB
dev.off()

figC <- ggplot(T1errorB %>% filter(Type == "Lorenz-96" & NZ == 0)) +
  geom_line(aes(x = Reference, y = Observed, color = Method, linetype = Method)) +
  scale_color_manual("Method", values = c("NPGC" = "#D62728", "cmlpGL" = "#1F77B4", "cmlpGSGL" = "#1F77B4",
                                          "cmlpH" = "#1F77B4", "RU" = "#FF7F0E", "GNS" = "#2CA02C")) +
  scale_linetype_manual("Method", values = c("NPGC" = "solid", "cmlpGL" = "dotted", "cmlpGSGL" = "dashed",
                                             "cmlpH" = "solid", "RU" = "solid", "GNS" = "solid")) +
  geom_abline(aes(slope = 1, intercept = 0), color = "lightgray") +
  facet_wrap(~ TimePts, nrow = 1,
             labeller = labeller(TimePts = c("250" = "T = 250", "500" = "T = 500", "1000" = "T = 1000"))) +
  scale_x_continuous("Level of Test", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  scale_y_continuous("Observed Type 1 Error", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank())

quartz(file = "/Users/noahgade/Desktop/Research/Rodu/NPGC/Figures/figureC.pdf", type = "pdf", width = 13, height = 5)
figC
dev.off()

figD <- ggplot(T1errorB %>% filter(Type == "Lorenz-96" & NZ == 2)) +
  geom_line(aes(x = Reference, y = Observed, color = Method, linetype = Method)) +
  scale_color_manual("Method", values = c("NPGC" = "#D62728", "cmlpGL" = "#1F77B4", "cmlpGSGL" = "#1F77B4",
                                          "cmlpH" = "#1F77B4", "RU" = "#FF7F0E", "GNS" = "#2CA02C")) +
  scale_linetype_manual("Method", values = c("NPGC" = "solid", "cmlpGL" = "dotted", "cmlpGSGL" = "dashed",
                                             "cmlpH" = "solid", "RU" = "solid", "GNS" = "solid")) +
  geom_abline(aes(slope = 1, intercept = 0), color = "lightgray") +
  facet_wrap(~ TimePts, nrow = 1,
             labeller = labeller(TimePts = c("250" = "T = 250", "500" = "T = 500", "1000" = "T = 1000"))) +
  scale_x_continuous("Level of Test", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  scale_y_continuous("Observed Type 1 Error", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank())

quartz(file = "/Users/noahgade/Desktop/Research/Rodu/NPGC/Figures/figureD.pdf", type = "pdf", width = 13, height = 5)
figD
dev.off()



# Table 3
GCident <- dplyr::tibble(Method = rep(c("NPGC", "cmlpGL", "cmlpGSGL", "cmlpH", "RU", "GNS"), each = 6), 
                         Type = rep(c("Lorenz-96", "TAR(2)"), 18),
                         TimePts = rep(rep(c(250, 500, 1000), each = 2), 6), GC = NA)

for(type in c("Lorenz-96", "TAR(2)")) {
  for(method in c("NPGC", "cmlpGL", "cmlpGSGL", "cmlpH", "RU", "GNS")) {
    temp250 <- sum(Results250 %>%
                     pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>% 
                     filter(Type == type & Causal == "GC" & Method == method) %>% select(Q) <= 0.1) 
    temp500 <- sum(Results500 %>%
                     pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>% 
                     filter(Type == type & Causal == "GC" & Method == method) %>% select(Q) <= 0.1) 
    temp1000 <- sum(Results1000 %>%
    pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>%
      filter(Type == type & Causal == "GC" & Method == method) %>% select(Q) <= 0.1) 
    
    GCident[(GCident$Type == type & GCident$Method == method & GCident$TimePts == 250),]$GC <-temp250
    GCident[(GCident$Type == type & GCident$Method == method & GCident$TimePts == 500),]$GC <-temp500
    GCident[(GCident$Type == type & GCident$Method == method & GCident$TimePts == 1000),]$GC <-temp1000
    
  }
}

GCident %>% mutate(GC = GC / 400) %>% pivot_wider(names_from = Method, values_from = GC)


GCidentB <- dplyr::tibble(Method = rep(c("NPGC", "cmlpGL", "cmlpGSGL", "cmlpH", "RU", "GNS"), each = 12), 
                         Type = rep(c("Lorenz-96", "TAR(2)"), 36),
                         TimePts = rep(rep(c(250, 500, 1000), each = 4), 6), 
                         NZ = rep(rep(c(0, 2), each = 2), 18), GC = NA)

for(type in c("Lorenz-96", "TAR(2)")) {
  for(method in c("NPGC", "cmlpGL", "cmlpGSGL", "cmlpH", "RU", "GNS")) {
    for(nz in c(0,2)) {
      temp250 <- sum(Results250 %>%
                       pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>% 
                       filter(Type == type & Causal == "GC" & Method == method & NZ == nz) %>% 
                       select(Q) <= 0.1) 
      temp500 <- sum(Results500 %>%
                       pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>% 
                       filter(Type == type & Causal == "GC" & Method == method& NZ == nz) %>% 
                       select(Q) <= 0.1) 
      temp1000 <- sum(Results1000 %>%
                        pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>%
                        filter(Type == type & Causal == "GC" & Method == method& NZ == nz) %>% 
                        select(Q) <= 0.1) 
      
      GCidentB[(GCidentB$Type == type & GCidentB$Method == method & GCidentB$TimePts == 250 & GCidentB$NZ == nz),]$GC <- temp250
      GCidentB[(GCidentB$Type == type & GCidentB$Method == method & GCidentB$TimePts == 500 & GCidentB$NZ == nz),]$GC <- temp500
      GCidentB[(GCidentB$Type == type & GCidentB$Method == method & GCidentB$TimePts == 1000 & GCidentB$NZ == nz),]$GC <- temp1000
    }
  }
}

GCidentB %>% mutate(GC = GC / 200) %>% pivot_wider(names_from = Method, values_from = GC)


load("~/Desktop/Research/Rodu/NPGC/Simulation Results/naive1000b.RData")
load("~/Desktop/Research/Rodu/NPGC/Simulation Results/npgc1000b.RData")
GL1000b <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/GL1000b.csv")
GSGL1000b <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/GSGL1000b.csv")
H1000b <- read.csv("~/Desktop/Research/Rodu/NPGC/Simulation Results/H1000b.csv")

Results1000b <- dplyr::tibble(NPGC = unlist(lapply(npgc1000b, extractQnpgc)),
                             cmlpGL = extractCMLP(GL1000b),
                             cmlpGSGL = extractCMLP(GSGL1000b),
                             cmlpH = extractCMLP(H1000b),
                             RU = unlist(lapply(naive1000b, extractQres, len = 1000, n_initializations = 1000)),
                             GNS = unlist(lapply(naive1000b, extractQnoise, len = 1000, n_initializations = 1000)),
                             TimePts = 1000,
                             Type = "TAR(2)",
                             Causal = "NC",
                             NZ = c(rep(0, 200), rep(2, 200)))

T1berror <- dplyr::tibble(Reference = rep(seq(1/200, 199/200, 1/200), 12),
                         Method = rep(rep(c("NPGC", "cmlpGL", "cmlpGSGL", "cmlpH", "RU", "GNS"), each = 199), 2),
                         TimePts = 1000,
                         Type = "TAR(2)",
                         NZ = c(rep(0, 1194), rep(2, 1194)),
                         Observed = NA)

for(method in c("NPGC", "cmlpGL", "cmlpGSGL", "cmlpH", "RU", "GNS")) {
  for(nz in c(0,2)) {
    temp1000 <- Results1000b %>%
      pivot_longer(cols = 1:6, names_to = "Method", values_to = "Q") %>%
      filter(Method == method & NZ == nz) %>% select(Q)
    
    T1berror[(T1berror$Method == method & T1berror$NZ == nz), ]$Observed <- sapply(seq(1/200, 199/200, 1/200), function(x) sum(temp1000 <= x) / 200)
  }
}

fig3 <- ggplot(T1berror) +
  geom_line(aes(x = Reference, y = Observed, color = Method, linetype = Method)) +
  scale_color_manual("Method", values = c("NPGC" = "#D62728", "cmlpGL" = "#1F77B4", "cmlpGSGL" = "#1F77B4",
                                          "cmlpH" = "#1F77B4", "RU" = "#FF7F0E", "GNS" = "#2CA02C")) +
  scale_linetype_manual("Method", values = c("NPGC" = "solid", "cmlpGL" = "dotted", "cmlpGSGL" = "dashed",
                                             "cmlpH" = "solid", "RU" = "solid", "GNS" = "solid")) +
  facet_wrap(~ NZ, nrow = 1,
             labeller = labeller(NZ = c("0" = "NZ = 0", "2" = "NZ = 2"))) +
  geom_abline(aes(slope = 1, intercept = 0), color = "lightgray") +
  scale_x_continuous("Level of Test", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  scale_y_continuous("Observed Type 1 Error", expand = c(0.002, 0.002), limits = c(0, 1),
                     breaks = seq(0.1, 1, 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank())

quartz(file = "/Users/noahgade/Desktop/Research/Rodu/NPGC/Figures/figure3.pdf", type = "pdf", width = 9, height = 5)
fig3
dev.off()





