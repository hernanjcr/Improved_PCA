#First load required libraries and codes
require(MASS)
require(ggplot2)
require(ggbiplot)
require(ggforce)
require(reshape)

#Setting initial random state
set.seed(54321)

#Procedure following the example published in stackexchange:
#  
#  https://stats.stackexchange.com/questions/35029/construct-artificial-slightly-overlapping-data-for-pca-plot
#
#nVars = number of variables
#rot   = random rotation matrix
#n1    = number of samples in group 1
#n2    = number of samples in group 2
#eps   = Error SD should be small compared to the SDs
#x     = simulated data
#y     = rotated simulated data
nVars <- 5
rot <- qr.Q(qr(matrix(rnorm(nVars*nVars), nVars)))
sigma <- function(theta=0, lambda=c(1,1)) {
  cos.t <- cos(theta);
  sin.t <- sin(theta)
  a <- matrix(c(cos.t, sin.t, -sin.t, cos.t), ncol=2)
  t(a) %*% diag(lambda) %*% a
}

n1 <- 50
n2 <- 75
x <- rbind(mvrnorm(n1, c(-2,-1), sigma(0, c(1/2,1))),
           mvrnorm(n2, c(0,1), sigma(pi/3, c(1, 1/3))))
eps <- 0.25
x <- cbind(x, matrix(rnorm(dim(x)[1]*(nVars-2), sd=eps), ncol=nVars-2))
y <- x %*% rot

colnames(y) <- paste(rep("X", 5), 1:5, sep='')
row.names(y) <- c(paste(rep("S", n1), 1:n1, sep=""), paste(rep("R",n2), 1:n2, sep=""))
groups <- as.factor(c(rep("group1", n1), rep("group2",n2)))

#Plotting the original data set.
cat('Plotting the original data. The ellipse corresponde to the sd of multinormal random.')
xdf <- data.frame(x)
ggplot(data=xdf, aes(x=X1, y=X2, color=groups, shape=groups)) + geom_point() +
  theme(legend.direction="horizontal", legend.position="top") +
  geom_ellipse(aes(x0=-2,y0=-1,b=1/2,a=1,angle=0)) +
  geom_ellipse(aes(x0=0,y0=1,b=1,a=1/3,angle=pi/3)) + coord_fixed() +
  ggtitle('Original data')

#Procesing the newly simulated data using PCA
pcaSim <- prcomp(y, center = TRUE, scale. = TRUE)

cat('Ploting the data after PCA transformation\n')
pcaSim.g <- ggbiplot(pcaSim, obs.scale=1, var.scale=1, groups=groups, ellipse=TRUE,
                     var.axes = FALSE) + scale_color_discrete(name="")  +
  theme(legend.direction="horizontal", legend.position="top") +
  ggtitle('PCA transformed data')
print(pcaSim.g)

#Additing variables with linear combination of previously created data:
yLin <- y
cat('\nNew variable:\n\tX6 <- 1.0*X1 + 2.0*X2\n')
yLin <- cbind(yLin, 1.0*y[,1] + 2.0*y[,2])
cat('New variable:\n\tX7 <- -0.5*X3 + 0.25*X4\n')
yLin <- cbind(yLin, -0.5*y[,3] + 0.25*y[,4])
colnames(yLin) <- paste(rep("X", 7), 1:7, sep='')

#Processing, as before, this newly data frame.
pcaLin <- prcomp(yLin, center = TRUE, scale. = TRUE)

cat('\nEigenvalues for the initial simulated data :\n')
cat(pcaSim$sdev^2)
cat('\nAfter linear combination of some columns:\n')
cat(pcaLin$sdev^2)

pltdf <- data.frame(PC=c(1:5,1:7),y=c(pcaSim$sdev^2/max(pcaSim$sdev^2),
                                      pcaLin$sdev^2/max(pcaLin$sdev^2)),
                    clase=as.factor(c(rep('Original',5),rep('Linear vars.',7))))
pcaSim.sc <- ggplot(pltdf, aes(PC,y,color=clase)) + geom_line() + geom_point() +
  scale_color_discrete(name="")+ theme(legend.direction="horizontal",
                                       legend.position="top") + xlab('principal component number') +
  ylab('proportion of explained variance') + scale_shape_manual("", values=c(19,15))
print(pcaSim.sc)

#Searching by columns with linear combinations as some eigenvalues are quasi-zero.
#There are two approach: First, follow the book "Principal Component Analisys" by Jolliffe and "Applied Multivariate Statistical Analisys" by K. Härdle and L. Simar.
cat('There are two small eigenvalues (6th and 7th).\n')
cat('The coefficients, loading values divided by the standard deviation, in the corresponding PC are:\n')
round(pcaLin$rotation[,6:7]/pcaLin$scale, digits = 2)
cat('The variable with the highest coefficient, in absolute value, is "X7"; as previously established X7 <- -0.5 * X3 + 0.25 * X4\n')
cat('From Härle and Simar book, calculating the correlation between variables and the respective PC :\n')
CorYpY <- pcaLin$rotation[,6:7]*pcaLin$sdev[6:7]/diag(cor(yLin))
round(apply(CorYpY, 2, function(x){return (x/max(x))}), digits=2)
cat('')

#Removing the identified columns
cat('Removing that variable and recalculing the PCA:\n')
pcaLin <- prcomp(yLin[,-7], center = TRUE, scale. = TRUE)
cat('\nThe new eigenvalues are :\n')
cat(pcaLin$sdev^2)
pcaSim.sc <- pcaSim.sc + geom_point(data =
                                      data.frame(x2=1:6,y2=pcaLin$sdev^2/max(pcaLin$sdev^2),clase=rep('without X7',6)),
                                    aes(x=x2,y=y2,colour=clase))
print(pcaSim.sc)
cat('\nThe coefficients for the 6th PC are now:\n')
round(pcaLin$rotation[,6]/pcaLin$scale, digits = 2)
cat('The Correlation of the variables to this PC is\n')
CorYpY2 <- pcaLin$rotation[,6]*pcaLin$sdev[6]/diag(cor(yLin[,-7]))
round(CorYpY2/max(abs(CorYpY2)), digits=2)

cat('On this case, the linear relaship can be clearly observed:\n
remembering that X6 is "X6 <- 1.0*X1 + 2.0*X2"\nor 0 = -X6 + X1 + 2*X2\n
equal to: X1 + 2*X2 + 0*X3 + 0*X4 + 0*X5 -X6\n
equal to the coefficients of the 6th PC, dividided by 0.24!')

cat('Eleminating the 2th variable, variable with the greatest coefficient, and recalculating the PCA:')
pcaLin <- prcomp(yLin[,-c(2,7)], center = TRUE, scale. = TRUE)
cat('\nThe new eigenvalues are now:\n')
cat(pcaLin$sdev^2)
cat('\nNow the lowest eigenvalue is 5% of the highest, that is, there is no linear relationship between the remaining variables.')
pcaSim.sc <- pcaSim.sc + geom_point(data =
                                      data.frame(x2=1:5,y2=pcaLin$sdev^2/max(pcaLin$sdev^2),clase=rep('Without X7 and X2',5)),
                                    aes(x=x2,y=y2, colour=clase))
print(pcaSim.sc)

cat('Ploting the final PCA.\n')
pcaSim.final.g <- ggbiplot(pcaLin, obs.scale=1, var.scale=1, groups=groups,
                           ellipse=TRUE, var.axes = FALSE) +
  scale_color_discrete(name="")  +
  theme(legend.direction="horizontal",legend.position="top")+
  ggtitle('PCA with linear relationships removed')
print(pcaSim.final.g)

