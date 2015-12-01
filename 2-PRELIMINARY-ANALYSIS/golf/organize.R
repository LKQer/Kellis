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

# Test different datasets: Cell Types
# train_ind <- (dat$bl$basic_chromosome == '18')
test_ind <- (dat$bl$basic_celltype == 'K562')
train_ind <- !test_ind
newdat <- list()
newy <- list()
newdat$train$all <- dat$all[train_ind,]
newdat$test$all <- dat$all[test_ind,]
newy$train$bin <- y$bin[train_ind]
newy$test$bin <- y$bin[test_ind]
newy$train$all <- y$all[train_ind,]

# Loading
load('dat.RData')
load('y.RData')

# Random Forest
library(randomForest)
rf <- randomForest(newdat$train$all, y = newy$train$bin, xtest = newdat$test$all, ytest = newy$test$bin, do.trace = 1, ntree = 5, keep.forest = TRUE)
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
pdf('echo.hist.pdf')
plot(p1, col=rgb(0, 0, 1, 1/2))
plot(p2, col=rgb(1, 0, 0, 1/2), add=T)
legend("topleft", c("Foreground", "Background"), col=c("blue", "red"), lwd=10)
dev.off()

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
