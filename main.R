
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
legend(58900, 50500, legend=c("ALL", "AML"),
       col=c("red", "blue"), lty=1:2, cex=0.8)

# prediction plot
#plot(p1, ncomp = nd, asp = 1, line = TRUE, type="n")
#text(Y_train, p1$fitted.values[,,nd], labels=rownames(X_train), col=as.vector(factor(train_response,levels=c(0,1),labels=c("red","blue"))))


predict_leukemia <- function(model, X_pls, threshold=0.5, ncomp=NULL) {
  if (!is.null(ncomp)) {
    predicted <- predict(model, newdata=X_pls, ncomp=ncomp, type="response")
  } else {
    predicted <- predict(model, newdata=X_pls, type="response")
  }
  
  predicted <- sapply(as.numeric(predicted), function(i) {return(if(i <= threshold) 0 else 1)})
  return(predicted)
}

thresholds <- seq(0,1, by=0.1)
accuracyTrain <- rep(0,length(thresholds))
accuracyTest <- rep(0,length(thresholds))
for (i in 1:length(thresholds)) {
  t <- thresholds[i]
  
  train_predicted <- predict_leukemia(p1, X_train, threshold=t, ncomp=nd)
  accuracyTrain[i] <- sum(train_predicted == train_response)/length(train_response)*100
  
  test_predicted <- predict_leukemia(p1, X_test, threshold=t, ncomp=nd)
  accuracyTest[i] <- sum(test_predicted == test_response)/length(test_predicted)*100
}

plot(rep(thresholds,2), c(accuracyTrain,accuracyTest),
     main="Accuracy with PSLR model", xlab="Thresholds", ylab="Accuracy")
lines(thresholds, accuracyTrain, col="orange", lwd=4)
lines(thresholds, accuracyTest, col="green", lwd=4)

# 4

X_train_pls <- p1$scores[,1:nd]
train_pls <- as.data.frame(cbind(X_train_pls, train_response))
X_train_pls <- train_pls[,which(colnames(train_pls)!="train_response")]

X_test_centered <- scale(X_test, center = colMeans(X_train), 
                    scale = FALSE)

X_test_pls <- X_test_centered %*% p1$projection
X_test_pls <- as.data.frame(X_test_pls[,1:nd])


# 5
plot(p1, plottype = "scores", comps = 1:2, type="n", main="X Scores")
text(X_train_pls, labels=rownames(X_train_pls), col=as.vector(factor(train_response,levels=c(0,1),labels=c("red","blue"))))
text(X_test_pls, labels=rownames(X_test_pls), col=as.vector(factor(test_response,levels=c(0,1),labels=c("darkred","darkblue"))))
legend(53100, 50500, legend=c("ALL train", "AML train", "ALL test", "AML test"),
       col=c("red", "blue", "darkred", "darkblue"), lty=c(1,1,1,1), cex=0.8)


# 6

logit.mod <- glm(train_response ~ ., data=train_pls, family=binomial(link="logit"))


# 7

thresholds <- seq(0,1, by=0.1)
accuracyTrain <- rep(0,length(thresholds))
accuracyTest <- rep(0,length(thresholds))
for (i in 1:length(thresholds)) {
  t <- thresholds[i]
  
  train_predicted <- predict_leukemia(logit.mod, X_train_pls, threshold=t)
  accuracyTrain[i] <- sum(train_predicted == train_response)/length(train_response)*100
  
  test_predicted <- predict_leukemia(logit.mod, X_test_pls, threshold=t)
  accuracyTest[i] <- sum(test_predicted == test_response)/length(test_predicted)*100
}

plot(rep(thresholds,2), c(accuracyTrain,accuracyTest),
     main="Accuracy with Logistic model", xlab="Thresholds", ylab="Accuracy")
lines(thresholds, accuracyTrain, col="orange", lwd=4)
lines(thresholds, accuracyTest, col="green", lwd=4)




# plot logistic predictions

summary(logit.mod)

df <- data.frame(c=X_train_pls[,1])
df$resp <- predict(logit.mod, newdata=X_train_pls, type="response")+0.00

df2 <- data.frame(c=X_test_pls[,1])
df2$resp <- predict(logit.mod, newdata=X_test_pls, type="response")-0.00

plot(df$c,df$resp, ylab="Response", xlab="Component 1", ylim=c(-0.05, 1.05), pch=16, cex=0.8,
     col=as.vector(factor(train_response,levels=c(0,1),labels=c("red","blue"))))

range <- range(df$c)
axe.x = seq(range[1],range[2],length=1000)
f.x = exp(fit$coef[1]+axe.x*fit$coef[2])/(1+exp(fit$coef[1]+axe.x*fit$coef[2]))
lines(axe.x,f.x,col="green",lwd=2)

points(df$c,df$resp, ylab="Response", xlab="Component 1", pch=16, cex=0.8,
       col=as.vector(factor(train_response,levels=c(0,1),labels=c("red","blue"))))
points(df2$c,df2$resp, ylab="Response", xlab="Component 1", pch=16, cex=0.8,
       col=as.vector(factor(test_response,levels=c(0,1),labels=c("darkred","darkblue"))))

legend(43500, 0.9, legend=c("ALL train", "AML train", "ALL test", "AML test", "link function"),
       col=c("red", "blue", "darkred", "darkblue", "green"), pch=c(16,16,16,16,-1), lty=c(0,0,0,0,1), cex=0.8)


