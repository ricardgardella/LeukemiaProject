
library("pls")

# 1 and 2

train <- read.csv("data/data_set_ALL_AML_train.csv", sep=";", header=TRUE)
test  <- read.csv("data/data_set_ALL_AML_independent.csv", sep=";", header=TRUE)

filter_dataset <- function(df) {
  numeric.colums <- gtools::mixedsort(colnames(df)[grepl("X", colnames(df))])
  df <- df[,numeric.colums]
  dft <- t(df)
  colnames(dft) <- paste0("g", seq(1,nrow(df)))
  dft <- as.data.frame(dft)
  return(dft)
}

traint <- filter_dataset(train)
testt <- filter_dataset(test)

train_response <- c(rep(0,27),rep(1,11))
test_response <- c(rep(0,11),rep(1,5),rep(0,2),rep(1,2),rep(0,1),rep(1,7),rep(0,6))

traint <- cbind(traint, train_response)
ind.resp <- which(colnames(traint)=="train_response")
X_train <- traint[,-ind.resp]
Y_train <- traint[,ind.resp]

X_test <- testt[,-ind.resp]


# 3

p1 <- plsr(train_response ~ ., center = TRUE, ncomp = 10, data = traint, validation = "LOO")
plot(RMSEP(p1), legendpos = "topright")

R2(p1)
plot(R2(p1), legendpos = "bottomright")

nd <- 4

plot(p1, plottype = "scores", comps = 1:2, type="n", main="X Scores")
text(p1$scores, labels=rownames(p1$scores), col=as.vector(factor(train_response,levels=c(0,1),labels=c("red","blue"))))
abline(h=0,v=0, col="gray")

# prediction plot
#plot(p1, ncomp = nd, asp = 1, line = TRUE, type="n")
#text(Y_train, p1$fitted.values[,,nd], labels=rownames(X_train), col=as.vector(factor(train_response,levels=c(0,1),labels=c("red","blue"))))


# 4

X_train_pls <- p1$scores
train_pls <- as.data.frame(cbind(X_train_pls, train_response))
X_train_pls <- train_pls[,which(colnames(train_pls)!="train_response")]

X_test_centered <- scale(X_test, center = colMeans(X_train), 
                    scale = FALSE)

X_test_pls <- X_test_centered %*% p1$projection
X_test_pls <- as.data.frame(X_test_pls)


# 5
plot(p1, plottype = "scores", comps = 1:2, type="n", main="X Scores")
text(X_train_pls, labels=rownames(X_train_pls), col=as.vector(factor(train_response,levels=c(0,1),labels=c("red","blue"))))
text(X_test_pls, labels=rownames(X_test_pls), col=as.vector(factor(test_response,levels=c(0,1),labels=c("darkred","darkblue"))))



# 6

logit.mod <- glm(train_response ~ ., data=train_pls, family=binomial(link="logit"))

predict_leukemia <- function(model, X_pls, threshold=0.5) {
  predicted <- predict(model, X_pls, type="response")
  predicted <- sapply(as.numeric(predicted), function(i) {return(if(i <= threshold) 0 else 1)})
  return(predicted)
}
train_predicted <- predict_leukemia(logit.mod, X_train_pls, threshold=0.7)
accuracy <- sum(train_predicted == train_response)/length(train_response)
accuracy*100


# 7

test_predicted <- predict_leukemia(logit.mod, X_test_pls, threshold=0.7)
accuracy = sum(test_predicted == test_response)/length(test_predicted)
accuracy*100




