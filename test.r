library(reticulate)
use_python("/home/aghee/anaconda3/bin/python")
np <- import("numpy")

cc <- c(94.91442, 39.03943, 75.73977, 59.46458, 62.45311, 58.24749)

jpeg("rplot.jpg", width = 350, height = 350)
plot(2:7, cc,
ylab="core consistency",
xlab="Number of components",
type="l",
ylim=c(0,100))
dev.off()

np$save("test.npy", cc)