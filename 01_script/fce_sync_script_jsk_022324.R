getwd()
sync <- read.csv("00_data/Smith_Fluorescence_Means_2019_2022_022324.csv")
sync
attach(sync)
View(sync)
head(sync)

library(synchrony)
#https://cran.r-project.org/web/packages/synchrony/synchrony.pdf

###Example from 'synchrony' package###
data(pisco.data)
View(pisco.data)
d=subset(pisco.data, subset=year==2000, select=c("latitude", "longitude", "sst"))
semiv=vario(data=d)
moran=vario(data=d, type="moran", nrands=100)
par(mfrow=c(2,1), mar=c(4.2, 4, 1, 1))
plot(semiv$mean.bin.dist, semiv$vario, xlab="Lag distance (km)", ylab="Semivariance")
plot(moran$mean.bin.dist, moran$vario, xlab="Lag distance (km)", ylab="Moran's I", t="l")
points(moran$mean.bin.dist[moran$pvals >= 0.05], moran$vario[moran$pvals >= 0.05],
       bg="white", pch=21)
points(moran$mean.bin.dist[moran$pvals < 0.05], moran$vario[moran$pvals < 0.05],
       bg="black", pch=21)
abline(h=0, lty=2)

# Compute spatial synchrony
d.upw=subset(pisco.data, select=c("latitude", "longitude", "year", "upwelling"))
d.cov=subset(pisco.data, select=c("latitude", "longitude", "year", "mussel_abund"))
# Reshape the data
d.upw.wide=reshape(data=d.upw, idvar=c("latitude", "longitude"), timevar=c("year"),
                   direction="wide")
d.cov.wide=reshape(data=d.cov, idvar=c("latitude", "longitude"), timevar=c("year"),
                   direction="wide")
# Generate variograms
v.upw=vario(n.bins=12, data=d.upw.wide, type="pearson", extent=1, nrands=999)

v.cov=vario(n.bins=12, data=d.cov.wide, type="pearson", extent=1, nrands=999)
## Fit variograms
v.cov.per=vario.fit(v.cov$vario, v.cov$mean.bin.dist, type="period",
                    start.vals=list(a=1, b=3, c=0))
v.upw.lin=vario.fit(v.upw$vario, v.upw$mean.bin.dist, type="linear")
par(mfrow=c(2,1))
plot(v.cov, xlab="Lag distance (km)", bg.sig="red", col.nonsig="red",
     main="Mussel cover",
     rug=TRUE, ylim=c(-0.3, 0.3))
lines(v.cov$mean.bin.dist, v.cov.per$fit, col="red")
plot(v.upw, xlab="Lag distance (km)", bg.sig="blue", col.nonsig="blue",
     main="Upwelling", rug=TRUE)
lines(v.upw$mean.bin.dist, v.upw.lin$fit, col="blue")
################################


#SRS vs TS/Ph Phase Synchrony
sync.srs=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C1", "Percent_C2", "Percent_C3", "Percent_C4", "Percent_C5", "Percent_C6")) 
sync.ts=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C1", "Percent_C2", "Percent_C3", "Percent_C4", "Percent_C5", "Percent_C6")) 


#C1-6 PARAFAC Components ALL
C1.sync=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C1")) 
semiv=vario(data=C1.sync, type="pearson")
moran=vario(data=C1.sync, type="pearson", nrands=1000)
plot(semiv$mean.bin.dist, semiv$vario, xlab="Lag distance (km)", ylab="Semivariance")
plot(moran$mean.bin.dist, moran$vario, xlab="Lag distance (km)", ylab="Moran's I", t="l")
points(moran$mean.bin.dist[moran$pvals >= 0.05], moran$vario[moran$pvals >= 0.05],
       bg="white", pch=21)
points(moran$mean.bin.dist[moran$pvals < 0.05], moran$vario[moran$pvals < 0.05],
       bg="black", pch=21)
abline(h=0, lty=2)

C2.sync=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C2")) 
semiv.C2.sync=vario(data=C2.sync,type="pearson")
moran.C2.sync=vario(data=C2.sync, type="pearson", nrands=1000)
par(mfrow=c(2,1), mar=c(4.2, 4, 1, 1))
plot(semiv.C2.sync$mean.bin.dist, semiv.C2.sync$vario,xlab="Lag distance (km)", ylab="Pearson") 
points(semiv.C2.sync$mean.bin.dist[moran.C2.sync$pvals < 0.05], moran.C2.sync$vario[moran.C2.sync$pvals < 0.05],
       bg="black", pch=21)
#abline(v = mean(semiv.sync), col="black", h=0, lty=2)

C3.sync=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C3")) 
semiv.C3.sync=vario(data=C3.sync,type="pearson")
moran.C3.sync=vario(data=C3.sync, type="pearson", nrands=1000)
par(mfrow=c(2,1), mar=c(4.2, 4, 1, 1))
plot(semiv.C3.sync$mean.bin.dist, semiv.C3.sync$vario,xlab="Lag distance (km)", ylab="Pearson") 
points(semiv.C3.sync$mean.bin.dist[moran.C3.sync$pvals < 0.05], moran.C3.sync$vario[moran.C3.sync$pvals < 0.05],
       bg="black", pch=21)
#abline(v = mean(semiv.sync), col="black", h=0, lty=2)

C4.sync=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C4")) 
semiv.C4.sync=vario(data=C4.sync,type="pearson")
moran.C4.sync=vario(data=C4.sync, type="pearson", nrands=1000)
par(mfrow=c(2,1), mar=c(4.2, 4, 1, 1))
plot(semiv.C4.sync$mean.bin.dist, semiv.C4.sync$vario,xlab="Lag distance (km)", ylab="Pearson") 
points(semiv.C4.sync$mean.bin.dist[moran.C4.sync$pvals < 0.05], moran.C4.sync$vario[moran.C4.sync$pvals < 0.05],
       bg="black", pch=21)
#abline(v = mean(semiv.sync), col="black", h=0, lty=2)

C5.sync=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C5")) 
semiv.C5.sync=vario(data=C5.sync,type="pearson")
moran.C5.sync=vario(data=C5.sync, type="pearson", nrands=1000)
par(mfrow=c(2,1), mar=c(4.2, 4, 1, 1))
plot(semiv.C5.sync$mean.bin.dist, semiv.C5.sync$vario,xlab="Lag distance (km)", ylab="Pearson") 
points(semiv.C5.sync$mean.bin.dist[moran.C5.sync$pvals < 0.05], moran.C5.sync$vario[moran.C5.sync$pvals < 0.05],
       bg="black", pch=21)
#abline(v = mean(semiv.sync), col="black", h=0, lty=2)

C6.sync=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C6")) 
semiv.C6.sync=vario(data=C6.sync,type="pearson")
moran.C6.sync=vario(data=C6.sync, type="pearson", nrands=1000)
par(mfrow=c(2,1), mar=c(4.2, 4, 1, 1))
plot(semiv.C6.sync$mean.bin.dist, semiv.C6.sync$vario,xlab="Lag distance (km)", ylab="Pearson") 
points(semiv.C6.sync$mean.bin.dist[moran.C6.sync$pvals < 0.05], moran.C6.sync$vario[moran.C6.sync$pvals < 0.05],
       bg="black", pch=21)
#abline(v = mean(semiv.sync), col="black", h=0, lty=2)



# Compute spatial synchrony
C1=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C1"))
C2=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C2"))
C3=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C3"))
C4=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C4"))
C5=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C5"))
C6=subset(sync, select=c("Year", "Latitude", "Longitude", "Percent_C6"))
str(C1)
# Reshape the data
C1.wide=reshape(data=C1, idvar=c("Latitude", "Longitude"), timevar=c("Year"), direction="wide")
C2.wide=reshape(data=C2, idvar=c("Latitude", "Longitude"), timevar=c("Year"), direction="wide")
C3.wide=reshape(data=C3, idvar=c("Latitude", "Longitude"), timevar=c("Year"), direction="wide")
C4.wide=reshape(data=C4, idvar=c("Latitude", "Longitude"), timevar=c("Year"), direction="wide")
C5.wide=reshape(data=C5, idvar=c("Latitude", "Longitude"), timevar=c("Year"), direction="wide")
C6.wide=reshape(data=C6, idvar=c("Latitude", "Longitude"), timevar=c("Year"), direction="wide")

# Generate variograms
v.C1=vario(n.bins=12, data=C1.wide, type="pearson", extent=1, nrands=9999)
v.C2=vario(n.bins=12, data=C2.wide, type="pearson", extent=1, nrands=9999)
v.C3=vario(n.bins=12, data=C3.wide, type="pearson", extent=1, nrands=9999)
v.C4=vario(n.bins=12, data=C4.wide, type="pearson", extent=1, nrands=9999)
v.C5=vario(n.bins=12, data=C5.wide, type="pearson", extent=1, nrands=9999)
v.C6=vario(n.bins=12, data=C6.wide, type="pearson", extent=1, nrands=9999)

## Fit variograms
v.C1.lin=vario.fit(v.C1$vario, v.C1$mean.bin.dist, type="linear")
v.C2.lin=vario.fit(v.C2$vario, v.C1$mean.bin.dist, type="linear")
v.C3.lin=vario.fit(v.C3$vario, v.C1$mean.bin.dist, type="linear")
v.C4.lin=vario.fit(v.C4$vario, v.C1$mean.bin.dist, type="linear")
v.C5.lin=vario.fit(v.C5$vario, v.C1$mean.bin.dist, type="linear")
v.C6.lin=vario.fit(v.C6$vario, v.C1$mean.bin.dist, type="linear")

par(mar=c(5,5,3,3))
par(mfrow=c(2,3))
plot(v.C1,  bg.sig="grey", col.nonsig="grey", xlab="Lag distance (km)",
     main="", rug=TRUE, ylim=c(-1, 1), cex.main=2, cex=2, cex.axis=2, cex.lab=2)
plot(v.C2,  bg.sig="grey", col.nonsig="grey", xlab="Lag distance (km)",
     main="", rug=TRUE, ylim=c(-1, 1), cex.main=2, cex=2, cex.axis=2, cex.lab=2)
plot(v.C3,  bg.sig="grey", col.nonsig="grey", xlab="Lag distance (km)",
     main="", rug=TRUE, ylim=c(-1, 1), cex.main=2, cex=2, cex.axis=2, cex.lab=2)
plot(v.C4,  bg.sig="grey", col.nonsig="grey", xlab="Lag distance (km)",
     main="", rug=TRUE, ylim=c(-1, 1), cex.main=2, cex=2, cex.axis=2, cex.lab=2)
plot(v.C5,  bg.sig="grey", col.nonsig="grey", xlab="Lag distance (km)",
     main="", rug=TRUE, ylim=c(-1, 1), cex.main=2, cex=2, cex.axis=2, cex.lab=2)
plot(v.C6,  bg.sig="grey", col.nonsig="grey", xlab="Lag distance (km)",
     main="", rug=TRUE, ylim=c(-1, 1), cex.main=2, cex=2, cex.axis=2, cex.lab=2)

