# Gets 72% test accuracy on toy dataset with 1000 interactions, 500 strong/500 weak
# Simple linear regression model

data <- read.table('/broad/compbio/maxwshen/data/1-MAKETRAINTEST/combined/toy.txt', header = TRUE)
y <- read.table('/broad/compbio/maxwshen/data/1-MAKETRAINTEST/combined/toy.y.txt', header = TRUE)

data <- cbind(data, y)

smp_size <- 750

set.seed(1)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]

model <- lm(Labels ~ ., train)
p <- predict(model, test)
acc <- sum((p > 1) == (test["Labels"] > 1))
total <- dim(test)[1]
cat("Test Accuracy:", acc / total, '\n')