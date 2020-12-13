
coord_mat <- as.matrix(expand.grid(x=1:9, y=1:9))
fmat <- matrix(rnorm(nrow(coord_mat)*100), nrow(coord_mat), 100)
cmat <- cor(t(fmat))
sfmat <- t(scale(t(fmat)))
wsa1 <- weighted_spatial_adjacency(coord_mat, fmat, nnk=3, weight_mode="binary", alpha=1, stochastic=TRUE)
wsa2 <- weighted_spatial_adjacency(coord_mat, sfmat, nnk=81, wsigma=.3,weight_mode="heat", alpha=0, normalized=FALSE,
                                   stochastic=FALSE, dthresh=50, sigma=25)


cor(as.matrix(wsa2)[lower.tri(wsa2)], cmat[lower.tri(cmat)])
plot(as.matrix(wsa2)[lower.tri(wsa2)], cmat[lower.tri(cmat)])
