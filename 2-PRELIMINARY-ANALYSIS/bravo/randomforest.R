library(clusterSim)
library(randomForest)

dat.imr90 <- read.table('~/Kellis/data/bravo/data.IMR90.txt', header = TRUE)
dat.gm12878 <- read.table('~/Kellis/data/bravo/data.GM12878.txt', header = TRUE)
dat.k562 <- read.table('~/Kellis/data/bravo/data.K562.txt', header = TRUE)

y.imr90 <- read.table('~/Kellis/data/bravo/y.IMR90.txt', header = TRUE)
y.gm12878 <- read.table('~/Kellis/data/bravo/y.GM12878.txt', header = TRUE)
y.k562 <- read.table('~/Kellis/data/bravo/y.K562.txt', header = TRUE)

dat.all <- merge(dat.imr90, dat.gm12878, all = TRUE)
dat.all <- merge(dat.all, dat.k562, all = TRUE)
dat.all <- dat.all[,colSums(is.na(dat.all)) == 0]
dim(dat.all)

y.all <- rbind(y.imr90, y.gm12878)
y.all <- rbind(y.all, y.k562)
y.bin <- (y.all > 1)
dim(y.all)

dists <- data.frame(dat.all[,"basic_genomic_dist"])

# Remove uninformative features
blacklist <- c("basic_chromosome", "basic_celltype")
dat.blacklist <- dat.all[, which(names(dat.all) %in% blacklist)]
dat.all <- dat.all[ , -which(names(dat.all) %in% blacklist)]

# Remove 0-variance columns from the data
temp <- lapply(dat.all, function(x) length(unique(x)))
novar <- which(!temp > 1)
dat.all <- dat.all[, - novar]

# Normalize to mean = 0, unit variance
dat.all <- data.Normalization(dat.all, type="n1", normalization="column")

# Split into trainig/test 
tr_size <- dim(dat.all)[1] * 0.80
set.seed(1)
train_ind <- sample(seq_len(nrow(dat.all)), size = tr_size)

dat.train <- dat.all[train_ind, ]
dat.test <- dat.all[-train_ind, ]
y.train_bin <- factor(y.bin[train_ind, ])
y.test_bin <- factor(y.bin[-train_ind, ])


# Build Random Forest model
rf <- randomForest(dat.train, y = y.train_bin, xtest = dat.test, ytest = y.test_bin, do.trace = 1, ntree = 151, keep.forest = TRUE)
rfp <- (rf$test$predicted > 1)
cat('Error Rate', sum(y.test_bin != rfp))
save(rf, file = 'rf.iter151.RData')