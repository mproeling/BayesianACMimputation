rm(list = ls(all = TRUE))

library(spatialprobit)

data("Katrina")
attach(Katrina) 
table(y1) # 300 of the 673 firms reopened during 0-3 months horizon, p.1016 
table(y2) # 425 of the 673 firms reopened during 0-6 months horizon, p.1016 
table(y3) # 478 of the 673 firms reopened during 0-12 months horizon, p.1016 detach(Katrina)
detach(Katrina)

if (require(ggmap)) { qmplot(long, lat, data = Katrina, maptype="roadmap", source="google") }
require(spdep)


# (a) 0-3 months time horizon # LeSage et al. (2011) use k=11 nearest neighbors in this case 
nb <- knn2nb(knearneigh(cbind(Katrina$lat, Katrina$long), k=11)) 
listw <- nb2listw(nb, style="W")
W1 <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")


fit1 <- sarprobit(y1 ~ flood_depth + log_medinc + small_size + large_size + low_status_customers + 
                    high_status_customers + owntype_sole_proprietor + owntype_national_chain, 
                  W=W1, data=Katrina, ndraw=600, burn.in = 100, showProgress=TRUE)
summary(fit1)

# (b) 0-6 months time horizon # LeSage et al. (2011) use k=15 nearest neighbors 
nb <- knn2nb(knearneigh(cbind(Katrina$lat, Katrina$long), k=15)) 
listw <- nb2listw(nb, style="W") 
W2 <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")

fit2 <- sarprobit(y2 ~ flood_depth + log_medinc + small_size + large_size + low_status_customers + 
                    high_status_customers + owntype_sole_proprietor + owntype_national_chain, 
                  W=W2, data=Katrina, ndraw=600, burn.in = 100, showProgress=TRUE) 
summary(fit2)

# LeSage et al. (2011) use k=15 nearest neighbors as in 0-6 months 
W3 <- W2
fit3 <- sarprobit(y3 ~ flood_depth + log_medinc + small_size + large_size + low_status_customers + 
                    high_status_customers + owntype_sole_proprietor + owntype_national_chain, 
                  W=W3, data=Katrina, ndraw=600, burn.in = 100, showProgress=TRUE) 
summary(fit3)

# replicate LeSage et al. (2011), Table 4, p.1018 # SAR probit model effects estimates for the 0-3-month time horizon 
impacts(fit1)
# replicate LeSage et al. (2011), Table 5, p.1019 # SAR probit model effects estimates for the 0-6-month time horizon 
impacts(fit2)
# replicate LeSage et al. (2011), Table 6, p.1020 # SAR probit model effects estimates for the 0-12-month time horizon 
impacts(fit3)

# example to build a spatial stochastic matrix using k nearest neighbour
W <- kNearestNeighbors(x=rnorm(100), y=rnorm(100), k=6)
image(W1, main="spatial weight matrix W")




