#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages # It's a library, so shhh!
shhh(library(tidyverse))
shhh(library(magrittr))

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied.", call. = FALSE)
}
limit <- as.numeric(args[1])

gaus <- suppressMessages(read_csv("gaussian_sampling.dat", col_names = c("Randoms")))
x <- 0.05
part <- pnorm(limit, mean = 0, sd = 1, lower.tail = FALSE)
gaus_true <- tibble(xs = seq(-5, 5, by = 0.01)) %>%
    mutate(ys = dnorm(xs))

plt <- gaus %>% ggplot(aes(x = Randoms, y = part*..density..)) + 
    geom_histogram(binwidth = x) +
    geom_line(data = gaus_true, aes(x = xs, y = ys), colour = "red") + 
    labs(x = "Value", y = "PDF")

ggsave("Gussian_tail_Sampling.png", plt)
