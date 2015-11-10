# Gets 81.2% test accuracy on alpha
# SVM Regression

library(e1071)  # cRAN SVM library

data_dir <- '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/combined/alpha/'
data_dir <- '~/Kellis/data/alpha/'    # ec2
data <- read.table(paste(data_dir, 'data.txt', sep = ''), header = TRUE)
y <- read.table(paste(data_dir, 'y.txt', sep = ''), header = TRUE)

data <- cbind(data, y)

smp_size <- 750
set.seed(1)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]

model <- svm(Labels ~ ., train)
p <- predict(model, test)
acc <- sum((p > 1) == (test["Labels"] > 1))
total <- dim(test)[1]
cat("Test Accuracy:", acc / total, '\n')


# Test rank of matrix
# qr(matrix)$rank