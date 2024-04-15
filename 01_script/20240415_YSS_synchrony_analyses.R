
# packages ----------------------------------------------------------------


if(!require(tidyverse)) install.packages("tidyverse", dependencies = TRUE)
if(!require(synchrony)) install.packages("synchrony", dependencies = TRUE)

############## data manipulation - ABUNDANCE -------------------------------------------------------


biota_coor <- readr::read_csv("00_data/biota_all_variables-2024-04-04.csv") |> 
  select(Location, Plot, Treatment, Longitude, Latitude) |> 
  rename(Site = Location) |> 
  mutate(Plot = as.numeric(Plot)) |> 
  unique()
biota_coor

names(biota_coor)

biota <- readr::read_csv("00_data/biota_seedlings.csv") |> 
  #filter(Species == "Euterpe edulis") |> 
  pivot_longer(cols = 5:24) |> 
  mutate(name = as.numeric(str_extract(name, "\\b\\d{4}\\b"))) |> 
  select(c(1:3,10,11)) |> 
  #filter(name == "T0" | name == "T108") |> 
  group_by(Site, Plot, Treatment, name) |> 
  summarise(sum = sum(value), .groups = "drop") |> 
  ungroup() |> 
  left_join(biota_coor) |> 
  rename(Time = name)  |> 
  pivot_wider(names_from = Time, values_from = sum) |> 
  select(!15) |> 
  as.data.frame()
biota



############## ABUNDANCE - Total ###############

# open ------------------------------------------------------------------


biota.open <- biota |> 
  filter(Treatment == "open") |> 
  # mutate(Treatment = str_replace(Treatment, "open", "0"),
  #        Treatment = str_replace(Treatment, "open", "1"),
  #        Treatment = as.numeric(Treatment)) |> 
  pivot_longer(cols = 6:14) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  #select(!c(1:3)) |> 
  as.data.frame()

str(biota.open)

biota.open

# Compute spatial synchrony
open <- subset(biota.open, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
open

# Reshape the data
open.wide <- reshape(data=open, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
open.wide

# Generate variograms
v.open <- vario(n.bins=10, data=open.wide, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.open

## Fit variograms
v.lin.open <- vario.fit(v.open$vario, v.open$mean.bin.dist, type="linear")
v.lin.open

# closed ------------------------------------------------------------------

biota.cld <- biota |> 
  filter(Treatment == "closed") |> 
  # mutate(Treatment = str_replace(Treatment, "closed", "0"),
  #        Treatment = str_replace(Treatment, "closed", "1"),
  #        Treatment = as.numeric(Treatment)) |> 
  pivot_longer(cols = 6:14) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  #select(!c(1:3)) |> 
  as.data.frame()

str(biota.cld)

biota.cld

# Compute spatial synchrony
closed <- subset(biota.cld, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
closed

# Reshape the data
closed.wide <- reshape(data=closed, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
closed.wide

# Generate variograms
v.closed <- vario(n.bins=10, data=closed.wide, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.closed

## Fit variograms
v.lin.closed <- vario.fit(v.closed$vario, v.closed$mean.bin.dist, type="linear")
v.lin.closed

# graphs ------------------------------------------------------------------

par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.open,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
#lines(v.open$mean.bin.dist, v.car.lin.open$fit, col="red")
plot(v.closed,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
mtext("Abundance", side = 3, cex.main=1.2, cex=1.2, line = - 2, outer = TRUE)



par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.lin.open,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
#lines(v.open$mean.bin.dist, v.car.lin.open$fit, col="red")
plot(v.lin.closed,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
mtext("Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)



############## ABUNDANCE - CBO ###############

# open ------------------------------------------------------------------


biota.cld.cbo <- biota |> 
  filter(Site == "CBO") |>
  filter(Treatment == "open") |> 
  # mutate(Treatment = str_replace(Treatment, "open", "0"),
  #        Treatment = str_replace(Treatment, "open", "1"),
  #        Treatment = as.numeric(Treatment)) |> 
  pivot_longer(cols = 6:14) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  #select(!c(1:3)) |> 
  as.data.frame()

str(biota.cld.cbo)

biota.cld.cbo

# Compute spatial synchrony
open.cbo <- subset(biota.cld.cbo, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
open.cbo

# Reshape the data
open.wide.cbo <- reshape(data=open.cbo, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
open.wide.cbo

# Generate variograms
v.open.cbo <- vario(n.bins=10, data=open.wide.cbo, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.open.cbo

## Fit variograms
v.cbo.lin.open.cbo <- vario.fit(v.open.cbo$vario, v.open.cbo$mean.bin.dist, type="linear")
v.cbo.lin.open.cbo

# closed ------------------------------------------------------------------

biota.cld.cbo <- biota |> 
  filter(Site == "CBO") |>
  filter(Treatment == "closed") |> 
  # mutate(Treatment = str_replace(Treatment, "closed", "0"),
  #        Treatment = str_replace(Treatment, "closed", "1"),
  #        Treatment = as.numeric(Treatment)) |> 
  pivot_longer(cols = 6:14) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  #select(!c(1:3)) |> 
  as.data.frame()

str(biota.cld.cbo)

biota.cld.cbo

# Compute spatial synchrony
closed.cbo <- subset(biota.cld.cbo, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
closed.cbo

# Reshape the data
closed.wide.cbo <- reshape(data=closed.cbo, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
closed.wide.cbo

# Generate variograms
v.closed.cbo <- vario(n.bins=10, data=closed.wide.cbo, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.closed.cbo

## Fit variograms
v.cbo.lin.closed.cbo <- vario.fit(v.closed.cbo$vario, v.closed.cbo$mean.bin.dist, type="linear")
v.cbo.lin.closed.cbo

# graphs ------------------------------------------------------------------

par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.open.cbo,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
#lines(v.open$mean.bin.dist, v.cbo.lin.open$fit, col="red")
plot(v.closed.cbo,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
mtext("CBO - Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)



par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.cbo.lin.open.cbo,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
#lines(v.open$mean.bin.dist, v.cbo.lin.open$fit, col="red")
plot(v.cbo.lin.closed.cbo,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
mtext("CBO - Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)



############## ABUNDANCE - VGM ###############

# open ------------------------------------------------------------------


biota.cld.vgm <- biota |> 
  filter(Site == "VGM") |>
  filter(Treatment == "open") |> 
  # mutate(Treatment = str_replace(Treatment, "open", "0"),
  #        Treatment = str_replace(Treatment, "open", "1"),
  #        Treatment = as.numeric(Treatment)) |> 
  pivot_longer(cols = 6:14) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  #select(!c(1:3)) |> 
  as.data.frame()

str(biota.cld.vgm)

biota.cld.vgm

# Compute spatial synchrony
open.vgm <- subset(biota.cld.vgm, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
open.vgm

# Reshape the data
open.wide.vgm <- reshape(data=open.vgm, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
open.wide.vgm

# Generate variograms
v.open.vgm <- vario(n.bins=10, data=open.wide.vgm, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.open.vgm

## Fit variograms
v.vgm.lin.open.vgm <- vario.fit(v.open.vgm$vario, v.open.vgm$mean.bin.dist, type="linear")


# closed ------------------------------------------------------------------

biota.cld.vgm <- biota |> 
  filter(Site == "VGM") |>
  filter(Treatment == "closed") |> 
  # mutate(Treatment = str_replace(Treatment, "closed", "0"),
  #        Treatment = str_replace(Treatment, "closed", "1"),
  #        Treatment = as.numeric(Treatment)) |> 
  pivot_longer(cols = 6:14) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  #select(!c(1:3)) |> 
  as.data.frame()

str(biota.cld.vgm)

biota.cld.vgm

# Compute spatial synchrony
closed.vgm <- subset(biota.cld.vgm, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
closed.vgm

# Reshape the data
closed.wide.vgm <- reshape(data=closed.vgm, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
closed.wide.vgm

# Generate variograms
v.closed.vgm <- vario(n.bins=10, data=closed.wide.vgm, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.closed.vgm

## Fit variograms
v.vgm.lin.closed.vgm <- vario.fit(v.closed.vgm$vario, v.closed.vgm$mean.bin.dist, type="linear")
v.vgm.lin.closed.vgm

# graphs ------------------------------------------------------------------

par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.open.vgm,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
#lines(v.open$mean.bin.dist, v.vgm.lin.open$fit, col="red")
plot(v.closed.vgm,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
mtext("VGM - Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)



par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.vgm.lin.open.vgm,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
#lines(v.open$mean.bin.dist, v.vgm.lin.open$fit, col="red")
plot(v.vgm.lin.closed.vgm,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
mtext("VGM - Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)



############## ABUNDANCE - ITA ###############

# open ------------------------------------------------------------------


biota.cld.ita <- biota |> 
  filter(Site == "ITA") |>
  filter(Treatment == "open") |> 
  # mutate(Treatment = str_replace(Treatment, "open", "0"),
  #        Treatment = str_replace(Treatment, "open", "1"),
  #        Treatment = as.numeric(Treatment)) |> 
  pivot_longer(cols = 6:14) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  #select(!c(1:3)) |> 
  as.data.frame()

str(biota.cld.ita)

biota.cld.ita

# Compute spatial synchrony
open.ita <- subset(biota.cld.ita, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
open.ita

# Reshape the data
open.wide.ita <- reshape(data=open.ita, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
open.wide.ita

# Generate variograms
v.open.ita <- vario(n.bins=10, data=open.wide.ita, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.open.ita

## Fit variograms
v.ita.lin.open.ita <- vario.fit(v.open.ita$vario, v.open.ita$mean.bin.dist, type="linear")
v.ita.lin.open.ita

# closed ------------------------------------------------------------------

biota.cld.ita <- biota |> 
  filter(Site == "ITA") |>
  filter(Treatment == "closed") |> 
  # mutate(Treatment = str_replace(Treatment, "closed", "0"),
  #        Treatment = str_replace(Treatment, "closed", "1"),
  #        Treatment = as.numeric(Treatment)) |> 
  pivot_longer(cols = 6:14) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  #select(!c(1:3)) |> 
  as.data.frame()

str(biota.cld.ita)

biota.cld.ita

# Compute spatial synchrony
closed.ita <- subset(biota.cld.ita, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
closed.ita

# Reshape the data
closed.wide.ita <- reshape(data=closed.ita, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
closed.wide.ita

# Generate variograms
v.closed.ita <- vario(n.bins=10, data=closed.wide.ita, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.closed.ita

## Fit variograms
v.ita.lin.closed.ita <- vario.fit(v.closed.ita$vario, v.closed.ita$mean.bin.dist, type="linear")
v.ita.lin.closed.ita

# graphs ------------------------------------------------------------------

par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.open.ita,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
#lines(v.open$mean.bin.dist, v.ita.lin.open$fit, col="red")
plot(v.closed.ita,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
mtext("ITA - Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)



par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.ita.lin.open.ita,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
#lines(v.open$mean.bin.dist, v.ita.lin.open$fit, col="red")
plot(v.ita.lin.closed.ita,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
mtext("ITA - Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)



############## ABUNDANCE - CAR ###############

# open ------------------------------------------------------------------


biota.cld.car <- biota |> 
  filter(Site == "CAR") |>
  filter(Treatment == "open") |> 
  # mutate(Treatment = str_replace(Treatment, "open", "0"),
  #        Treatment = str_replace(Treatment, "open", "1"),
  #        Treatment = as.numeric(Treatment)) |> 
  pivot_longer(cols = 6:14) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  #select(!c(1:3)) |> 
  as.data.frame()

str(biota.cld.car)

biota.cld.car

# Compute spatial synchrony
open.car <- subset(biota.cld.car, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
open.car

# Reshape the data
open.wide.car <- reshape(data=open.car, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
open.wide.car

# Generate variograms
v.open.car <- vario(n.bins=10, data=open.wide.car, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.open.car

## Fit variograms
v.car.lin.open.car <- vario.fit(v.open.car$vario, v.open.car$mean.bin.dist, type="linear")
v.car.lin.open.car

# closed ------------------------------------------------------------------

biota.cld.car <- biota |> 
  filter(Site == "CAR") |>
  filter(Treatment == "closed") |> 
  # mutate(Treatment = str_replace(Treatment, "closed", "0"),
  #        Treatment = str_replace(Treatment, "closed", "1"),
  #        Treatment = as.numeric(Treatment)) |> 
  pivot_longer(cols = 6:14) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  #select(!c(1:3)) |> 
  as.data.frame()

str(biota.cld.car)

biota.cld.car

# Compute spatial synchrony
closed.car <- subset(biota.cld.car, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
closed.car

# Reshape the data
closed.wide.car <- reshape(data=closed.car, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
closed.wide.car

# Generate variograms
v.closed.car <- vario(n.bins=10, data=closed.wide.car, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.closed.car

## Fit variograms
v.car.lin.closed.car <- vario.fit(v.closed.car$vario, v.closed.car$mean.bin.dist, type="linear")
v.car.lin.closed.car

# graphs ------------------------------------------------------------------

par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.open.car,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
#lines(v.open$mean.bin.dist, v.car.lin.open$fit, col="red")
plot(v.closed.car,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
mtext("CAR - Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)



par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.car.lin.open.car,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
#lines(v.open$mean.bin.dist, v.car.lin.open$fit, col="red")
plot(v.car.lin.closed.car,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
mtext("CAR - Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)



############## data manipulation - DIVERSITY  -------------------------------------------------------


biota_coor <- readr::read_csv("00_data/biota_all_variables-2024-04-04.csv") |> 
  select(Location, Plot, Treatment, Longitude, Latitude) |> 
  rename(Site = Location) |> 
  mutate(Plot = as.numeric(Plot)) |> 
  unique()
biota_coor

names(biota_coor)

biota.div <- readr::read_csv("00_data/Simp_div_area.csv") |> 
  #filter(Species == "Euterpe edulis") |> 
  pivot_longer(cols = 4:13) |> 
  mutate(name = as.numeric(str_extract(name, "\\b\\d{4}\\b"))) |> 
  pivot_wider(names_from = name, values_from = value) |> 
  left_join(biota_coor) |> 
  as.data.frame()
biota.div



############## DIVERSITY - Total ###############
# open ------------------------------------------------------------------


biota.open.div <- biota.div |> 
  filter(Treatment == "open") |> 
  pivot_longer(cols = 4:13) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  as.data.frame()

str(biota.open.div)

biota.open.div

# Compute spatial synchrony
open.div <- subset(biota.open.div, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
open.div

# Reshape the data
open.wide.div <- reshape(data=open.div, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
open.wide.div

# Generate variograms
v.open.div <- vario(n.bins=10, data=open.wide.div, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.open.div

## Fit variograms
v.lin.open.div <- vario.fit(v.open.div$vario, v.open.div$mean.bin.dist, type="linear")
v.lin.open.div

# closed ------------------------------------------------------------------

biota.cld.div <- biota.div |> 
  filter(Treatment == "closed") |> 
  pivot_longer(cols = 4:13) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  as.data.frame()

str(biota.cld.div)

biota.cld.div

# Compute spatial synchrony
closed.div <- subset(biota.cld.div, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
closed.div

# Reshape the data
closed.wide.div <- reshape(data=closed.div, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
closed.wide.div

# Generate variograms
v.closed.div <- vario(n.bins=10, data=closed.wide.div, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.closed.div

## Fit variograms
v.lin.closed.div <- vario.fit(v.closed.div$vario, v.closed.div$mean.bin.dist, type="linear")
v.lin.closed.div

# graphs ------------------------------------------------------------------

par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.open.div,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
#lines(v.open$mean.bin.dist, v.car.lin.open$fit, col="red")
plot(v.closed.div,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
mtext("Simpson Diversity - Species", side = 3, cex.main=1.2, cex=1.2, line = - 2, outer = TRUE)



par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.lin.open.div,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
#lines(v.open$mean.bin.dist, v.car.lin.open$fit, col="red")
plot(v.lin.closed.div,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
mtext("Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)



############## data manipulation - DIVERSITY  -------------------------------------------------------


biota_coor <- readr::read_csv("00_data/biota_all_variables-2024-04-04.csv") |> 
  select(Location, Plot, Treatment, Longitude, Latitude) |> 
  rename(Site = Location) |> 
  mutate(Plot = as.numeric(Plot)) |> 
  unique()
biota_coor

names(biota_coor)

biota.div.lf <- readr::read_csv("00_data/invs_simp_among.csv") |> 
  pivot_longer(cols = 4:13) |> 
  mutate(name = as.numeric(str_extract(name, "\\b\\d{4}\\b"))) |> 
  pivot_wider(names_from = name, values_from = value) |> 
  left_join(biota_coor) |> 
  as.data.frame()
biota.div.lf



############## DIVERSITY Life forms - Total ###############
# open ------------------------------------------------------------------


biota.open.div.lf <- biota.div.lf |> 
  filter(Treatment == "open") |> 
  pivot_longer(cols = 4:13) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  as.data.frame()

str(biota.open.div.lf)

biota.open.div.lf

# Compute spatial synchrony
open.div.lf <- subset(biota.open.div.lf, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
open.div.lf

# Reshape the data
open.wide.div.lf <- reshape(data=open.div.lf, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
open.wide.div.lf

# Generate variograms
v.open.div.lf <- vario(n.bins=10, data=open.wide.div.lf, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.open.div.lf

## Fit variograms
v.lin.open.div.lf <- vario.fit(v.open.div.lf$vario, v.open.div.lf$mean.bin.dist, type="linear")
v.lin.open.div.lf

# closed ------------------------------------------------------------------

biota.cld.div.lf <- biota.div.lf |> 
  filter(Treatment == "closed") |> 
  pivot_longer(cols = 4:13) |> 
  rename(year = name) |> 
  mutate(year = as.numeric(year)) |> 
  as.data.frame()

str(biota.cld.div.lf)

biota.cld.div.lf

# Compute spatial synchrony
closed.div.lf <- subset(biota.cld.div.lf, select=c("year", "Latitude", "Longitude", "value")) |> drop_na()
closed.div.lf

# Reshape the data
closed.wide.div.lf <- reshape(data=closed.div.lf, idvar=c("Latitude", "Longitude"), timevar=c("year"), direction="wide")
closed.wide.div.lf

# Generate variograms
v.closed.div.lf <- vario(n.bins=10, data=closed.wide.div.lf, type="pearson", extent=1, nrands=500, size.bins=NULL)
v.closed.div.lf

## Fit variograms
v.lin.closed.div.lf <- vario.fit(v.closed.div.lf$vario, v.closed.div.lf$mean.bin.dist, type="linear")
v.lin.closed.div.lf

# graphs ------------------------------------------------------------------

par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.open.div.lf,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
#lines(v.open$mean.bin.dist, v.car.lin.open$fit, col="red")
plot(v.closed.div.lf,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1, ci = TRUE)
mtext("Simpson Diversity - Life-forms", side = 3, cex.main=1.2, cex=1.2, line = - 2, outer = TRUE)


par(mar=c(5,5,3,3))
par(mfrow=c(1,2))
plot(v.lin.open.div.lf,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Open", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
lines(v.open$mean.bin.dist, v.car.lin.open$fit, col="red")
plot(v.lin.closed.div.lf,  bg.sig="grey40", col.nonsig="grey20", xlab="Lag distance (km)",
     main="Closed", rug=TRUE, ylim=c(-1, 1), cex.main=1, cex=1, cex.axis=1, cex.lab=1,  ci = TRUE)
mtext("Abundance", side = 3, cex.main=1.5, cex=1.5, line = - 2, outer = TRUE)

