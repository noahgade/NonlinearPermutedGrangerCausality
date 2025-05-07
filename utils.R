library(tidyverse)
library(Matrix)
library(MASS)
library(cowplot)

npgc_eval <- function(npgc_result, alpha = 0.1) {
  # Return <= alpha, Granger Causal
  Output <- vector("list", length = 3)
  Output[[1]] <- mean(rowSums(npgc_result)[1] >= rowSums(npgc_result))
  Output[[2]] <- (Output[[1]] <= alpha)
  Output[[3]] <- alpha
  names(Output) <- c("Quantile", "Decision", "Level")
  return(Output)
}

gauss_eval <- function(gauss_result, alpha = 0.1) {
  # Return <= alpha, Granger Causal
  Output <- vector("list", length = 3)
  Output[[1]] <- mean(rowSums(gauss_result)[1] >= rowSums(gauss_result))
  Output[[2]] <- (Output[[1]][1] <= alpha)
  Output[[3]] <- alpha
  names(Output) <- c("Quantile", "Decision", "Level")
  return(Output)
}

zero_eval <- function(zero_result, prop = 0.5) {
  # Return >= prop, Granger Causal
  Output <- vector("list", length = 3)
  Output[[1]] <- mean(apply(zero_result, 2, function(zz) zz[1] <= zz[2]))
  Output[[2]] <- (Output[[1]][1] >= prop)
  Output[[3]] <- prop
  names(Output) <- c("Votes", "Decision", "Level")
  return(Output)
}

lasso_eval <- function(lasso_result, prop = 0.5, unique = FALSE) {
  # Return >= prop, Granger Causal
  Output <- vector("list", length = 4)
  Output[[1]] <- mean(colSums(lasso_result$W0) != 0)
  Output[[2]] <- (Output[[1]][1] >= prop)
  Output[[3]] <- prop
  Output[[4]] <- lasso_result$Lambdas
  if(unique) {
    Output[[4]] <- unique(Output[[4]])
  }
  names(Output) <- c("Votes", "Decision", "Level", "Lambda")
  return(Output)
}


# Functions to plot simulated data and neural network results

plotData <- function(data, times = NULL, vars = NULL, predict = NULL, error = FALSE, rescale = FALSE) {
  # Function to visualize data
  # data: output from data generating function including 
  #       $data, a matrix of data (L x Dim)
  #       $Cpts, change points in the data (optional)
  #       $States, states associated with the change points (optional)
  #       $Regimes, regimes associated with the state (optional)
  #       $Causal, a matrix of causal connections among variables (optional)
  # times: rows to plot, default is all
  # vars: columns to plot, default is all
  # predict: adds predicted values
  # error: if TRUE, plots the error (data - predict) instead of the data
  # rescale: if TRUE and error = TRUE, rescales the vertical axis to zoom in on the error
  L <- nrow(data$data)
  Dim <- ncol(data$data)
  if(is.null(times)) {
    times <- 1:L
  }
  if(is.null(vars)) {
    vars <- 1:Dim
  }
  p <- vector("list", length = length(vars))
  for(index in 1:length(vars)) {
    p[[index]] <- local({
      L <- L; Dim <- Dim; times <- times; vars <- vars; index <- index; 
      data <- data; predict <- predict; error <- error; rescale <- rescale
      start <- min(times); end <- max(times); xrange <- end - start + 1
      if(error) {
        to_plot <- data$data[times, vars, drop = FALSE] - predict
        if(rescale) {
          lower <- min(to_plot); upper <- max(to_plot)
        } else {
          lower <- min(data$data[times, vars, drop = FALSE]); upper <- max(data$data[times, vars, drop = FALSE])
        }
      } else {
        to_plot <- data$data[times, vars, drop = FALSE]
        lower <- min(to_plot); upper <- max(to_plot)
      }
      yrange = upper - lower
      cpts <- data$Cpts[,1][data$Cpts[,1] %in% times[1:(length(times) - 1)]]
      states <- data$States[times,]
      xticks <- 10 * ceiling(xrange / 100)
      plot <- ggplot() +
        scale_x_continuous(limits = c(start - 0.5, end + 0.5), expand = c(0, 0), breaks = seq(xticks, end - 1, xticks)) +
        scale_y_continuous(limits = c(lower - 0.05 * yrange, upper + 0.05 * yrange), expand = c(0, 0)) +
        scale_fill_manual(values = c("0" = alpha("gray100", 0.5), "1" = alpha("gray80", 0.5), "2" = alpha("gray60", 0.5), "R0" = "gray100", "R1" = "gray40"), guide = "none") +
        theme_minimal() +
        geom_rect(aes(xmin = c(start - 0.5, cpts + 0.5), 
                      xmax = c(cpts + 0.5, end + 0.5), 
                      ymin = rep(lower - 0.05 * yrange, length(cpts) + 1), 
                      ymax = rep(upper + 0.05 * yrange, length(cpts) + 1),
                      fill = factor(states[c(diff(states), 1) != 0]))) +
        geom_line(aes(x = start:end, y = to_plot[,index]), linewidth = 0.5)
      if(!is.null(predict) & !error) {
        plot <- plot + geom_line(aes(x = start:end, y = predict[,index]), linewidth = 0.8, color = "red", alpha = 0.5)
      }
      if(error & !rescale) {
        plot <- plot + annotate("segment", x = start, xend = end, y = lower, yend = lower, linewidth = 0.5, color = "gray40", linetype = "dotted") +
          annotate("segment", x = start, xend = end, y = upper, yend = upper, linewidth = 0.5, color = "gray40", linetype = "dotted")
      }
      if(!is.null(data$Regimes)) {
        regimes <- paste0("R", data$Regimes[times,])
        plot <- plot + geom_rect(aes(xmin = times - 0.5, xmax = times + 0.5, ymin = upper + 0.02 * yrange, ymax = upper + 0.05 * yrange, fill = regimes))
      }
      if((index == 1) & (length(vars) == 1))  {
        plot <- plot + theme(axis.text.x = element_text(size = 12),
                             axis.text.y = element_blank(),
                             axis.title = element_blank(),
                             axis.ticks = element_blank(),
                             plot.margin = unit(c(4, 4, 2, 2), "mm"),
                             panel.grid.major.y = element_blank(),
                             panel.grid.minor.y = element_blank(),
                             panel.grid.minor.x = element_blank(),
                             panel.border = element_rect(color = "white", fill = NA)) +
          geom_label(aes(x = cpts + 0.025 * xrange, y = lower + 0.15 * yrange, label = cpts), size = 5)
      } else if((index == 1) & (length(vars) > 1)) {
        plot <- plot + theme(axis.text = element_blank(),
                             axis.title = element_blank(),
                             axis.ticks = element_blank(),
                             plot.margin = unit(c(4, 4, -0.5, 2), "mm"),
                             panel.grid.major.y = element_blank(),
                             panel.grid.minor.y = element_blank(),
                             panel.grid.minor.x = element_blank(),
                             panel.border = element_rect(color = "white", fill = NA))
      } else if((index == length(vars)) & (length(vars) > 1))  {
        plot <- plot + theme(axis.text.x = element_text(size = 12),
                             axis.text.y = element_blank(),
                             axis.title = element_blank(),
                             axis.ticks = element_blank(),
                             plot.margin = unit(c(-0.5, 4, 2, 2), "mm"),
                             panel.grid.major.y = element_blank(),
                             panel.grid.minor.y = element_blank(),
                             panel.grid.minor.x = element_blank(),
                             panel.border = element_rect(color = "white", fill = NA)) +
          geom_label(aes(x = cpts + 0.025 * xrange, y = lower + 0.15 * yrange, label = cpts), size = 5)
      } else {
        plot <- plot + theme(axis.text = element_blank(),
                             axis.title = element_blank(),
                             axis.ticks = element_blank(),
                             plot.margin = unit(c(-0.5, 4, -0.5, 2), "mm"),
                             panel.grid.major.y = element_blank(),
                             panel.grid.minor.y = element_blank(),
                             panel.grid.minor.x = element_blank(),
                             panel.border = element_rect(color = "white", fill = NA))
      }
      print(plot)
    })
  }
  cowplot::plot_grid(plotlist = p, ncol = 1, rel_heights = c(rep(1, length(vars) - 1), 1.1))
}

produce_tableV1 <- function(alpha) {
  rows <- 1:12
  tablefull <- tibble(Length = rep(c(250, 500, 1000), each = 4), Model = rep(rep(c("VAR", "TAR"), each = 2), 3), Decision = rep(c("GC", "NC"), 6)) 
  tablefull <- tablefull %>% 
    mutate(NPGC25.2.tanh = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_1[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>% 
    mutate(NPGC25.2.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_2[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>% 
    mutate(NPGC25.2.leaky = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_3[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>% 
    mutate(NPGC25.4.tanh = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_4[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>% 
    mutate(NPGC25.4.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_5[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>% 
    mutate(NPGC25.4.leaky = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_6[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>% 
    mutate(NPGC25.8.tanh = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_7[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>% 
    mutate(NPGC25.8.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_8[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC25.8.leaky = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_9[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC50.2.tanh = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_10[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC50.2.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_11[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC50.2.leaky = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_12[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC50.4.tanh = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_13[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC50.4.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_14[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC50.4.leaky = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_15[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC50.8.tanh = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_16[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC50.8.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_17[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC50.8.leaky = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_18[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC100.2.tanh = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_19[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC100.2.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_20[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC100.2.leaky = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_21[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC100.4.tanh = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_22[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC100.4.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_23[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC100.4.leaky = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_24[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC100.8.tanh = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_25[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC100.8.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_26[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(NPGC100.8.leaky = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_27[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))}))
}

produce_tableV2 <- function(alpha) {
  rows <- 1:12
  tablefull <- tibble(Length = rep(c(250, 500, 1000), each = 4), Model = rep(rep(c("VAR", "TAR"), each = 2), 3), Decision = rep(c("GC", "NC"), 6)) 
  tablefull <- tablefull %>% 
    mutate(NPGC100.2.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_20[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(GAUSS = sapply(rows, function(zz) {mean(unlist(lapply(GAUSS100[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {gauss_eval(xx, alpha)$Decision})))})) %>% 
    mutate(ZERO = sapply(rows, function(zz) {mean(unlist(lapply(ZERO100[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {zero_eval(xx)$Decision})))}))
}

produce_tableV3 <- function(alpha) {
  rows <- 1:36
  tablefull <- tibble(Length = rep(c(250, 500, 1000), each = 12), Model = rep(rep(c("VAR", "TAR"), each = 6), 3), Decision = rep(c("GC1", "GC2", "NC1", "NC2", "NC3", "NC4"), 6)) 
  tablefull <- tablefull %>% 
    mutate(NPGC100.2.relu = sapply(rows, function(zz) {mean(unlist(lapply(NPGC_20[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {npgc_eval(xx, alpha)$Decision})))})) %>%
    mutate(GAUSS = sapply(rows, function(zz) {mean(unlist(lapply(GAUSS100[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {gauss_eval(xx, alpha)$Decision})))})) %>% 
    mutate(ZERO = sapply(rows, function(zz) {mean(unlist(lapply(ZERO100[which(str_detect(datakey$Model, tablefull$Model[zz]) & str_detect(datakey$GC, tablefull$Decision[zz]) & datakey$Length == tablefull$Length[zz])], function(xx) {zero_eval(xx)$Decision})))}))
}