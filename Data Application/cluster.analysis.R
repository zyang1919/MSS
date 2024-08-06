## clean history 
rm(list = ls())
library(MASS)
library(LaplacesDemon)
library(plotly)
library(dplyr)
library(plot3D)
library(Matrix)
library(car)
library(rgl)
library(cluster)

sur_dt <- "~/Desktop/Seemingly Unrelated Regression/Data"
dist.mall <- get(load(paste(sur_dt, "dist.alltime.all.(mean) 1.RData", sep = '/')))
dt <- get(load(paste(sur_dt, "MCMCestimation.mean 1.RData", sep = '/')))
tm <- c(0.0,  0.1,0.2,  0.3,  0.4,  0.5,  0.6,  0.7,0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9,2.0)



#############
## plot 
#############
#############
## cluster 1
#############
C1 <- get(load(paste(sur_dt, "Cluster=1.all.RData", sep = '/')))
c1.dt <- C1$data
c1.smooth <- C1$smooth
c1.time <- C1$time
scatter3D(x=c1.smooth[,1], y = c1.smooth[,2], z = c1.smooth[,3], colvar = xout,
          theta = 40, phi = 20, col = ramp.col(c("red", "green", "blue")),
          type='l', bty = "b2", lwd=2, xlab="Air Puff", ylab="Drop", zlab='Shake', lwd=2, 
          xlim = c(0, 3.5), ylim = c(0, 3.5), zlim = c(0,3.5), 
          ticktype = "detailed", clab = "Time")


############
## cluster 2
############
C2 <- get(load(paste(sur_dt, "Cluster=2.all.RData", sep = '/')))
c2.dt <- C2$data
c2.smooth <- C2$smooth
c2.time <- C2$time

lines3D(x=c2.smooth[,1], y = c2.smooth[,2], z = c2.smooth[,3], colvar = xout, 
        col = ramp.col(c("red", "green", "blue")), theta = 40, phi = 20, 
        type='l', bty = "b2", lwd=2, xlab="Air Puff", ylab="Drop", zlab='Shake', lwd=2, 
        xlim = c(0, 1.0), ylim = c(0.4, 1.0), zlim = c(0,1.0), 
        ticktype = "detailed", clab = "Time", add = TRUE)

####### t=0
d2.0 <- c2.time[[1]]
center2.0 <- c(mean(d2.0[,1]), mean(d2.0[,2]), mean(d2.0[,3]))
cov_matrix.2.0 <- cov(d2.0)
ellipsoid.2.0 <- ellipse3d(cov_matrix.2.0, centre = center2.0, scale = c(.1, .1, .1))
ellipsoid_x.2.0 <- ellipsoid.2.0$vb[1, ] / ellipsoid.2.0$vb[4, ]
ellipsoid_y.2.0 <- ellipsoid.2.0$vb[2, ] / ellipsoid.2.0$vb[4, ]
ellipsoid_z.2.0 <- ellipsoid.2.0$vb[3, ] / ellipsoid.2.0$vb[4, ]
scatter3D(ellipsoid_x.2.0, ellipsoid_y.2.0, ellipsoid_z.2.0, add = TRUE, col = "#FF0000FF", alpha = 1, pch = ".")

####### t=1
d2.1 <- c2.time[[2]]
center2.1 <- c(mean(d2.1[,1]), mean(d2.1[,2]), mean(d2.1[,3]))
cov_matrix.2.1 <- cov(d2.1)
ellipsoid.2.1 <- ellipse3d(cov_matrix.2.1, centre = center2.1, scale = c(.1, .1, .1))
ellipsoid_x.2.1 <- ellipsoid.2.1$vb[1, ] / ellipsoid.2.1$vb[4, ]
ellipsoid_y.2.1 <- ellipsoid.2.1$vb[2, ] / ellipsoid.2.1$vb[4, ]
ellipsoid_z.2.1 <- ellipsoid.2.1$vb[3, ] / ellipsoid.2.1$vb[4, ]
scatter3D(ellipsoid_x.2.1, ellipsoid_y.2.1, ellipsoid_z.2.1, add = TRUE, col = "#02FC00FF", alpha = 1, pch = ".")

####### t=2
d2.2 <- c2.time[[3]]
center2.2 <- c(mean(d2.2[,1]), mean(d2.2[,2]), mean(d2.2[,3]))
cov_matrix.2.2 <- cov(d2.2)
ellipsoid.2.2 <- ellipse3d(cov_matrix.2.2, centre = center2.2, scale = c(.1, .1, .1))
ellipsoid_x.2.2 <- ellipsoid.2.2$vb[1, ] / ellipsoid.2.2$vb[4, ]
ellipsoid_y.2.2 <- ellipsoid.2.2$vb[2, ] / ellipsoid.2.2$vb[4, ]
ellipsoid_z.2.2 <- ellipsoid.2.2$vb[3, ] / ellipsoid.2.2$vb[4, ]
scatter3D(ellipsoid_x.2.2, ellipsoid_y.2.2, ellipsoid_z.2.2, add = TRUE, col = "#0000FFFF", alpha = 1, pch = ".")
####################
### end of cluster 2 
####################


############
## cluster 3
############
C3 <- get(load(paste(sur_dt, "Cluster=3.all.RData", sep = '/')))
c3.dt <- C3$data
c3.smooth <- C3$smooth
c3.time <- C3$time

lines3D(x=c3.smooth[,1], y = c3.smooth[,2], z = c3.smooth[,3], colvar = xout, 
        col = ramp.col(c("red", "green", "blue")), theta = 40, phi = 20, 
        type='l', bty = "b2", lwd=2, xlab="Air Puff", ylab="Drop", zlab='Shake', lwd=2, 
        xlim = c(0, 1.0), ylim = c(0.4, 1.0), zlim = c(0,1.0), 
        ticktype = "detailed", clab = "Time", add = TRUE)

####### t=0
d3.0 <- c3.time[[1]]
center3.0 <- c(mean(d3.0[,1]), mean(d3.0[,2]), mean(d3.0[,3]))
cov_matrix.3.0 <- cov(d3.0)
ellipsoid.3.0 <- ellipse3d(cov_matrix.3.0, centre = center3.0, scale = c(.1, .1, .1))
ellipsoid_x.3.0 <- ellipsoid.3.0$vb[1, ] / ellipsoid.3.0$vb[4, ]
ellipsoid_y.3.0 <- ellipsoid.3.0$vb[2, ] / ellipsoid.3.0$vb[4, ]
ellipsoid_z.3.0 <- ellipsoid.3.0$vb[3, ] / ellipsoid.3.0$vb[4, ]
scatter3D(ellipsoid_x.3.0, ellipsoid_y.3.0, ellipsoid_z.3.0, add = TRUE, col = "#FF0000FF", alpha = 1, pch = ".")

####### t=1
d3.1 <- c3.time[[2]]
center3.1 <- c(mean(d3.1[,1]), mean(d3.1[,2]), mean(d3.1[,3]))
cov_matrix.3.1 <- cov(d3.1)
ellipsoid.3.1 <- ellipse3d(cov_matrix.3.1, centre = center3.1, scale = c(.1, .1, .1))
ellipsoid_x.3.1 <- ellipsoid.3.1$vb[1, ] / ellipsoid.3.1$vb[4, ]
ellipsoid_y.3.1 <- ellipsoid.3.1$vb[2, ] / ellipsoid.3.1$vb[4, ]
ellipsoid_z.3.1 <- ellipsoid.3.1$vb[3, ] / ellipsoid.3.1$vb[4, ]
scatter3D(ellipsoid_x.3.1, ellipsoid_y.3.1, ellipsoid_z.3.1, add = TRUE, col = "#02FC00FF", alpha = 1, pch = ".")

####### t=2
d3.2 <- c3.time[[3]]
center3.2 <- c(mean(d3.2[,1]), mean(d3.2[,2]), mean(d3.2[,3]))
cov_matrix.3.2 <- cov(d3.2)
ellipsoid.3.2 <- ellipse3d(cov_matrix.3.2, centre = center3.2, scale = c(.1, .1, .1))
ellipsoid_x.3.2 <- ellipsoid.3.2$vb[1, ] / ellipsoid.3.2$vb[4, ]
ellipsoid_y.3.2 <- ellipsoid.3.2$vb[2, ] / ellipsoid.3.2$vb[4, ]
ellipsoid_z.3.2 <- ellipsoid.3.2$vb[3, ] / ellipsoid.3.2$vb[4, ]
scatter3D(ellipsoid_x.3.2, ellipsoid_y.3.2, ellipsoid_z.3.2, add = TRUE, col = "#0000FFFF", alpha = 1, pch = ".")
####################
### end of cluster 3 
####################


############
## cluster 4
############
C4 <- get(load(paste(sur_dt, "Cluster=4.all.RData", sep = '/')))
c4.dt <- C4$data
c4.smooth <- C4$smooth
c4.time <- C4$time

lines3D(x=c4.smooth[,1], y = c4.smooth[,2], z = c4.smooth[,3], colvar = xout, 
        col = ramp.col(c("red", "green", "blue")), theta = 20, phi = 20, 
        type='l', bty = "b2", lwd=2, xlab="Air Puff", ylab="Drop", zlab='Shake', lwd=2, 
        xlim = c(0, 1.0), ylim = c(0.4, 1.0), zlim = c(0,1.0), 
        ticktype = "detailed", clab = "Time", add = TRUE)

####### t=0
d4.0 <- c4.time[[1]]
center4.0 <- c(mean(d4.0[,1]), mean(d4.0[,2]), mean(d4.0[,3]))
cov_matrix.4.0 <- cov(d4.0)
ellipsoid.4.0 <- ellipse3d(cov_matrix.4.0, centre = center4.0, scale = c(.1, .1, .1))
ellipsoid_x.4.0 <- ellipsoid.4.0$vb[1, ] / ellipsoid.4.0$vb[4, ]
ellipsoid_y.4.0 <- ellipsoid.4.0$vb[2, ] / ellipsoid.4.0$vb[4, ]
ellipsoid_z.4.0 <- ellipsoid.4.0$vb[3, ] / ellipsoid.4.0$vb[4, ]
scatter3D(ellipsoid_x.4.0, ellipsoid_y.4.0, ellipsoid_z.4.0, add = TRUE, col = "#FF0000FF", alpha = 1, pch = ".")

####### t=1
d4.1 <- c4.time[[2]]
center4.1 <- c(mean(d4.1[,1]), mean(d4.1[,2]), mean(d4.1[,3]))
cov_matrix.4.1 <- cov(d4.1)
ellipsoid.4.1 <- ellipse3d(cov_matrix.4.1, centre = center4.1, scale = c(.1, .1, .1))
ellipsoid_x.4.1 <- ellipsoid.4.1$vb[1, ] / ellipsoid.4.1$vb[4, ]
ellipsoid_y.4.1 <- ellipsoid.4.1$vb[2, ] / ellipsoid.4.1$vb[4, ]
ellipsoid_z.4.1 <- ellipsoid.4.1$vb[3, ] / ellipsoid.4.1$vb[4, ]
scatter3D(ellipsoid_x.4.1, ellipsoid_y.4.1, ellipsoid_z.4.1, add = TRUE, col = "#02FC00FF", alpha = 1, pch = ".")

####### t=2
d4.2 <- c4.time[[3]]
center4.2 <- c(mean(d4.2[,1]), mean(d4.2[,2]), mean(d4.2[,3]))
cov_matrix.4.2 <- cov(d4.2)
ellipsoid.4.2 <- ellipse3d(cov_matrix.4.2, centre = center4.2, scale = c(.1, .1, .1))
ellipsoid_x.4.2 <- ellipsoid.4.2$vb[1, ] / ellipsoid.4.2$vb[4, ]
ellipsoid_y.4.2 <- ellipsoid.4.2$vb[2, ] / ellipsoid.4.2$vb[4, ]
ellipsoid_z.4.2 <- ellipsoid.4.2$vb[3, ] / ellipsoid.4.2$vb[4, ]
scatter3D(ellipsoid_x.4.2, ellipsoid_y.4.2, ellipsoid_z.4.2, add = TRUE, col = "#0000FFFF", alpha = 1, pch = ".")
####################
### end of cluster 4 
####################