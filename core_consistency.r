library(multiway)
library(reticulate)
use_python("/home/aghee/anaconda3/bin/python")
np <- import("numpy")
X = np$load('results/sim1/mean.npy')
cmin = 2
cmax = 40
nstart = 20

cc <- matrix(0, cmax-cmin+1, nstart)
for (R in cmin:cmax){
    for (point in 1:nstart){
        pfac <- parafac(X,nfac=R,nstart=1, const = c("nonneg", "nonneg", "nonneg"))
        cc[R-cmin+1, point] <- corcondia(X, pfac)
    }
    print(cc)
}

# jpeg("results/core_consistency.jpg", width = 350, height = 350)
# plot(2:cmax, cc,
# ylab="core consistency",
# xlab="Number of components",
# type="l",
# ylim=c(0,100))
# dev.off()
np$save("results/sim1/ccmean.npy", cc)

