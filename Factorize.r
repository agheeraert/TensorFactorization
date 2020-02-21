library(multiway)
library(reticulate)
use_python("/home/aghee/anaconda3/bin/python")
np <- import("numpy")
X = np$load('results/apo_sim1.npy')
cmin = 3
cmax = 7

for (R in cmin:cmax){
    pfac <- parafac(X,nfac=R,nstart=20, const = c("nonneg", "nonneg", "nonneg"), output="best")
    np$save(paste0("results/ABC", R, ".npy"), pfac)
}

# jpeg("results/core_consistency.jpg", width = 350, height = 350)
# plot(2:cmax, cc,
# ylab="core consistency",
# xlab="Number of components",
# type="l",
# ylim=c(0,100))
# dev.off()

