#################################################
# Simulate data
#################################################

set.seed(2)
data <- data.frame(x = c(rnorm(600, 10, 2.7), rnorm(300, 20, 2.7)))
data <- data[sample.int(900),,drop=F]

library(ggplot2)
ggplot(data, aes(x = x)) + geom_histogram(binwidth = 0.5) + theme_bw() +
  theme(panel.border = element_blank(), axis.ticks = element_blank()) + ylab("Count") + xlab("Size") + ggtitle("Cat or Tiger?")

#################################################
# Build Model
#################################################

varDim <- c(2, 100)
facPot <- list(multinomialPotential(dim = c(1, 2)),
               normalPotential(dim = c(2, 100)))
facNbs <- list(c(1),
               c(1,2))
mixDfg <- dfg(varDim = varDim, facPot = facPot, facNbs = facNbs, varNames = c("H", "X"), facNames = c("Prior H", "X | H"))
plot(mixDfg)

#################################################
# Format data
#################################################

library(dplyr)
library(dplyr)
data_discrete <- data %>%
  mutate(H = NA) %>%
  mutate(X = as.integer(cut(x, breaks = seq(0,30,length.out = 101), labels = c(1:100)))) %>%
  select(H, X)

#################################################
# Train model
#################################################

optimFun <- list(mixNorm = normOptimize(range = c(0, 30)))
train(data_discrete, mixDfg, optim = c("row", "mixNorm"), optimFun = optimFun, iter.max = 500)

#################################################
# Inspect
#################################################

pot <- potentials(mixDfg)[[2]]
library(reshape2)
pot_df <- melt(pot)
pot_df$size <- seq(0, 30, length.out = 101)[pot_df$Var2]
pot_df$Var1 <- as.factor(pot_df$Var1)
ggplot(pot_df, aes(x = size, y = value, colour = as.factor(Var1))) + geom_line()

mixMps <- mps(data = data_discrete, dfg = mixDfg)
head(mixMps, 10)