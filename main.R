
setwd("C:\\Users\\ester\\Dropbox\\0. KMLMM\\Part 2\\Repo\\data")

train <- read.csv("data_set_ALL_AML_train.csv", sep=";", header=TRUE)
test  <- read.csv("data_set_ALL_AML_independent.csv", sep=";", header=TRUE)

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

response <- c(rep(1,27),rep(0,11))
traint <- cbind(traint, response)
ind.resp <- which(colnames(traint)=="response")
X_train <- traint[,-ind.resp]
Y_train <- traint[,ind.resp]

X_test <- test[,-ind.resp]

p1 <- plsr(response ~ ., ncomp = 10, data = traint, validation = "LOO")
plot(RMSEP(p1), legendpos = "topright")

R2(p1)
plot(R2(p1), legendpos = "bottomright")

nd <- 4

plot(p1, plottype = "scores", comps = 1:2, type="n", main="X Scores")
text(p1$scores, labels=rownames(p1$scores), col=as.vector(factor(response,levels=c(0,1),labels=c("red","blue"))))
abline(h=0,v=0, col="gray")

# prediction plot
plot(p1, ncomp = nd, asp = 1, line = TRUE, type="n")
text(Y_train, p1$fitted.values[,,nd], labels=rownames(X_train), col=as.vector(factor(response,levels=c(0,1),labels=c("red","blue"))))

