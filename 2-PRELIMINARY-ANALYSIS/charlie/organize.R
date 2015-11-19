library(clusterSim)
library(randomForest)

read_data <- function(typ, dy) {
  ct <- list()
  ct$imr90 <- read.table(paste('/broad/compbio/maxwshen/data/1-MAKETRAINTEST/combined/charlie/', typ, '/', dy, '.IMR90.txt', sep=''), header = TRUE)
  ct$k562 <- read.table(paste('/broad/compbio/maxwshen/data/1-MAKETRAINTEST/combined/charlie/', typ, '/', dy, '.K562.txt', sep=''), header = TRUE)
  return(ct)
}

merge_data <- function(dat) {
  m <- merge(dat$imr90, dat$k562, all = TRUE)
  m <- m[,colSums(is.na(m)) == 0]
  return(m)
}

merge_y <- function(y) {
  m <- list()
  m$all <- rbind(y$imr90, y$k562)
  m$bin <- factor(m$all > 1.5)
  return(m)
}

blacklist <- function(dat) {
  blacklist <- c("basic_chromosome", "basic_celltype")
  dat$bl <- dat$all[, which(names(dat$all) %in% blacklist)]
  dat$all <- dat$all[ , -which(names(dat$all) %in% blacklist)]
  return(dat)
}

normalize <- function(dat) {
  lentrain <- dim(dat$train$all)[1]
  m <- rbind(dat$train$all, dat$test$all)

  # Remove 0-variance columns from the data
  temp <- lapply(m, function(x) length(unique(x)))
  novar <- which(!temp > 1)
  m <- m[, - novar]

  # Mean-center the columns
  m <- as.data.frame(sapply(m, as.numeric))
  m <- scale(m)

  dat$train$all <- m[1: lentrain,]
  dat$test$all <- m[(lentrain + 1):dim(m)[1],]
  return(dat)
}

dat <- list()
y <- list()
dat$train <- read_data('train', 'data')
dat$test <- read_data('test', 'data')
y$train <- read_data('train', 'y')
y$test <- read_data('test', 'y')

dat$train$all <- merge_data(dat$train)
dat$test$all <- merge_data(dat$test)
y$train <- merge_y(y$train)
y$test <- merge_y(y$test)

dat$train <- blacklist(dat$train)
dat$test <- blacklist(dat$test)

dat <- normalize(dat)

save(dat, file = 'dat.RData')
save(y, file = 'y.RData')

# Build Random Forest model
rf <- randomForest(dat$train$all, y = y$train$bin, xtest = dat$test$all, ytest = y$test$bin, do.trace = 1, ntree = 45, keep.forest = TRUE)
rfp <- (rf$test$predicted > 1)
cat('Error Rate', sum(y.test_bin != rfp))
save(rf, file = 'rf.iter151.RData')