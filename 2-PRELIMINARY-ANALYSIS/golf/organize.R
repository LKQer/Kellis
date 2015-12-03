library(gtools)

read_data <- function(dy) {
  ct <- list()
  ct$imr90 <- read.table(paste('/broad/compbio/maxwshen/data/1-MAKETRAINTEST/combined/golf/', dy, '.IMR90.txt', sep=''), header = TRUE)
  ct$k562 <- read.table(paste('/broad/compbio/maxwshen/data/1-MAKETRAINTEST/combined/golf/', dy, '.K562.txt', sep=''), header = TRUE)
  ct$gm12878 <- read.table(paste('/broad/compbio/maxwshen/data/1-MAKETRAINTEST/combined/golf/', dy, '.GM12878.txt', sep=''), header = TRUE)
  return(ct)
}

merge_data <- function(dat) {
  # Merge different celltypes which have different columns
  # m <- merge(dat$imr90, dat$k562, all = TRUE)
  m <- smartbind(dat$imr90, dat$k562)
  m <- smartbind(m, dat$gm12878)
  m <- m[,colSums(is.na(m)) == 0]
  return(m)
}

merge_y <- function(y) {
  # Merge y in different celltypes
  m <- list()
  m$all <- rbind(y$imr90, y$k562)
  m$all <- rbind(m$all, y$gm12878)
  m$bin <- factor(m$all > 1.5)
  return(m)
}

blacklist <- function(dat) {
  blacklist <- c("basic_chromosome", "basic_celltype")
  dat$bl <- dat$all[, which(names(dat$all) %in% blacklist)]
  dat$all <- dat$all[ , -which(names(dat$all) %in% blacklist)]
  dat$bl$basic_genomic_dist <- dat$all$basic_genomic_dist
  return(dat)
}

normalize <- function(dat) {
  lentrain <- dim(dat$all)[1]
  m <- dat$all

  # Remove 0-variance columns from the data
  temp <- lapply(m, function(x) length(unique(x)))
  novar <- which(!temp > 1)
  m <- m[, - novar]

  # Mean-center the columns
  m <- as.data.frame(sapply(m, as.numeric))
  m <- scale(m)

  dat$all <- m
  return(dat)
}

# Data preparation Functions
# ------------------------------------------------------------------------------
# Data prep

dat <- list()
y <- list()
dat <- read_data('data')
y <- read_data('y')

dat$all <- merge_data(dat)
y <- merge_y(y)

dat <- blacklist(dat)
dat <- normalize(dat)

dat$imr90 <- dat$all[dat$bl$basic_celltype == 'IMR90',]
dat$k562 <- dat$all[dat$bl$basic_celltype == 'K562',]
dat$gm12878 <- dat$all[dat$bl$basic_celltype == 'GM12878',]
y$imr90 <- y$bin[dat$bl$basic_celltype == 'IMR90']
y$k562 <- y$bin[dat$bl$basic_celltype == 'K562']
y$gm12878 <- y$bin[dat$bl$basic_celltype == 'GM12878']

save(dat, file = 'dat.RData')
save(y, file = 'y.RData')

# Loading
load('dat.RData')
load('y.RData')

# Data preparation
# ------------------------------------------------------------------------------
# Analysis

# Test different datasets: Cell Types
distances <- (dat$bl$basic_genomic_dist > 0) & (dat$bl$basic_genomic_dist < 10000000)
# test_ind <- (dat$bl$basic_celltype == 'K562')
# train_ind <- (dat$bl$basic_celltype != 'K562' & dat$bl$basic_chromosome != '1')
# valid_ind <- (dat$bl$basic_celltype != 'K562' & dat$bl$basic_chromosome == '1')
test_ind <- (dat$bl$basic_chromosome == '1')
train_ind <- (dat$bl$basic_chromosome != '1' & dat$bl$basic_chromosome != '3')
valid_ind <- (dat$bl$basic_chromosome != '1' & dat$bl$basic_chromosome == '3')
test_ind <- test_ind & distances
train_ind <- train_ind & distances
newdat <- list()
newy <- list()
newdat$train$all <- dat$all[train_ind,]
newdat$test$all <- dat$all[test_ind,]
newdat$valid$all <- dat$all[valid_ind,]
newy$train$bin <- y$bin[train_ind]
newy$test$bin <- y$bin[test_ind]
newy$valid$bin <- y$bin[valid_ind]
newy$train$all <- y$all[train_ind,]

sum(test_ind); sum(train_ind); sum(valid_ind)

# Write.table to file
idd <- 'B'
write.table(newdat$train$all, file = paste(idd, '.train.dat.txt', sep = ''), sep='\t')
write.table(newy$train$bin, file = paste(idd, '.train.y.txt', sep = ''), sep='\t')
write.table(newdat$test$all, file = paste(idd, '.test.dat.txt', sep = ''), sep='\t')
write.table(newy$test$bin, file = paste(idd, '.test.y.txt', sep = ''), sep='\t')
write.table(newdat$valid$all, file = paste(idd, '.valid.dat.txt', sep = ''), sep='\t')
write.table(newy$valid$bin, file = paste(idd, '.valid.y.txt', sep = ''), sep='\t')


# Random Forest
library(randomForest)
rf <- randomForest(newdat$train$all, y = newy$train$bin, xtest = newdat$test$all, ytest = newy$test$bin, do.trace = 1, ntree = 21, keep.forest = TRUE)
save(rf, file = 'rf.iter151.RData')


# ROC
library(ROCR)
prob <- predict(rf, newdata = newdat$test$all, type='response')
pred <- prediction(as.numeric(prob), as.numeric(newy$test$bin))
perf <- performance(pred, measure = 'tpr', x.measure = 'fpr')

pdf('roc.pdf')
plot(perf, col=rainbow(10))
dev.off()

auc <- performance(pred, measure = 'auc')
auc <- auc@y.values[[1]]
cat('AUC:', auc, '\n')

# PCA
pc <- prcomp(dat$all)
library(ggbiplot)
# classes <- factor(y$bin, levels = c(FALSE, TRUE), labels = c('Background', 'Foreground'))
# classes <- dat$bl$basic_celltype
classes <- dat$bl$basic_genomic_dist
pdf('golf.pca.pdf')
# g <- ggbiplot(pc, obs.scale = 1, var.axes = FALSE, groups = classes, ellipse = TRUE, alpha = 0.1)
g <- ggbiplot(pc, obs.scale = 1, var.axes = FALSE, groups = classes, ellipse = FALSE, alpha = 0.1)
# g <- g + scale_color_discrete(name = '')
g <- g + scale_color_gradient(name = '')
print(g)
dev.off()

# Histogram
p1 <- hist(dat$bl$basic_genomic_dist[y$all > 1])
p2 <- hist(dat$bl$basic_genomic_dist[y$all < 1])
pdf('golf.hist.pdf')
plot(p1, col=rgb(0, 0, 1, 1/2))
plot(p2, col=rgb(1, 0, 0, 1/2), add=T)
legend("topleft", c("Foreground", "Background"), col=c("blue", "red"), lwd=10)
dev.off()

# Analysis
# ------------------------------------------------------------------------------
# ETC

# GLMNET
library(glmnet)
ytrainallfixed <- as.vector(sapply(y$train$all, as.numeric))
datfixed <- as.matrix(dat$train$all)
glm <- glmnet(datfixed, ytrainallfixed)
# do not know how to interpret...

# Linear Model
m <- lm(newy$train$all ~., data.frame(newdat$train$all))
p <- predict(m, data.frame(newdat$test$all))
sum((p > 1) != newy$test$bin)
