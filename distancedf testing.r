source("distancedf.r")
setwd("C:/Users/rowcliffe.m/OneDrive - Zoological Society of London/CameraTrapping/REM/Training/REM example/2017")
dat <- read.csv("BCI_firstdetection_data.csv")
summary(dat)
dat$humid[dat$humid==999] <- NA
dat$temp[dat$temp==999] <- NA
summary(dat)


m1 <- fitdf(radius~1, subset(dat, species=="agouti"), transect="point")
plot(m1$ddf, pdf=TRUE)
m2 <- fitdf(radius~1, subset(dat, species=="agouti"), transect="point", order=0)
plot(m2$ddf, pdf=TRUE)
m3 <- fitdf(radius~1, subset(dat, species=="agouti"), transect="point", key="hr")
plot(m3$ddf, pdf=TRUE)
m4 <- fitdf(radius~1, subset(dat, species=="agouti"), transect="point", key="unif")
plot(m4$ddf, pdf=TRUE)
m1$ddf$criterion
m3$ddf$criterion
m4$ddf$criterion
mods <- list(m1, m3, m4)
AICdf(mods)

m5 <- fitdf(radius~log(mass)+time+season, dat, transect="point", key="hr")
m5$edd
nd <- data.frame(mass=1, season="DRY", time="DAY")
m6 <- fitdf(radius~log(mass)+time+season, dat, transect="point", key="hr", newdata=nd)
m6$ddf$par
m6$edd
nd <- data.frame(mass=seq(0.1,30,len=10), season="DRY", time="NIGHT")
m7 <- fitdf(radius~log(mass)+time+season, dat, transect="point", key="hr", newdata=nd)
m7$edd

m8 <- fitdf(radius~log(mass)+time+season, dat, transect="point", key="hr")
m9 <- fitdf(radius~log(mass)+time, dat, transect="point", key="hr")
m10 <- fitdf(radius~log(mass)+season, dat, transect="point", key="hr")
m11 <- fitdf(radius~time+season, dat, transect="point", key="hr")
m12 <- fitdf(radius~log(mass), dat, transect="point", key="hr")
m13 <- fitdf(radius~time, dat, transect="point", key="hr")
m14 <- fitdf(radius~season, dat, transect="point", key="hr")
m15 <- fitdf(radius~1, dat, transect="point", key="hr", adjustment=NULL)
mods <- list(m8, m9, m10, m11, m12, m13, m14, m15)
AICdf(mods)

dat2 <- subset(dat, !is.na(temp) & !is.na(humid))
nd <- data.frame(expand.grid(mass=unique(dat2$mass), time=levels(dat2$time)), temp=mean(dat2$temp), humid=mean(dat2$humid))
m16 <- fitdf(radius~log(mass)+time+temp+humid, dat2, transect="point", key="hr", newdata=nd)
m17 <- fitdf(radius~log(mass)+time+temp, dat2, transect="point", key="hr", newdata=nd)
m18 <- fitdf(radius~log(mass)+time+humid, dat2, transect="point", key="hr", newdata=nd)
m19 <- fitdf(radius~log(mass)+temp+humid, dat2, transect="point", key="hr", newdata=nd)
m20 <- fitdf(radius~time+temp+humid, dat2, transect="point", key="hr", newdata=nd)
m21 <- fitdf(radius~log(mass)+time, dat2, transect="point", key="hr", newdata=nd)
m22 <- fitdf(radius~log(mass)+temp, dat2, transect="point", key="hr", newdata=nd)
m23 <- fitdf(radius~time+temp, dat2, transect="point", key="hr", newdata=nd)
m24 <- fitdf(radius~log(mass)+humid, dat2, transect="point", key="hr", newdata=nd)
m25 <- fitdf(radius~time+humid, dat2, transect="point", key="hr", newdata=nd)
m26 <- fitdf(radius~temp+humid, dat2, transect="point", key="hr", newdata=nd)
m27 <- fitdf(radius~log(mass), dat2, transect="point", key="hr", newdata=nd)
m28 <- fitdf(radius~time, dat2, transect="point", key="hr", newdata=nd)
m29 <- fitdf(radius~temp, dat2, transect="point", key="hr", newdata=nd)
m30 <- fitdf(radius~humid, dat2, transect="point", key="hr", newdata=nd)
m31 <- fitdf(radius~1, dat2, transect="point", key="hr")
mods2 <- list(m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31)
mm <- c(2, which(unlist(lapply(mods2, function(m) is.null(m$edd)))==TRUE))
AICdf(mods2[-mm])

#Plotting responses
i <- order(unique(dat$mass))
pdat <- subset(m6$edd, season=="LATEWET" & time=="DAY")
cl=2
with(pdat, plot(mass[i], estimate[i], type="l", log="xy", col=cl, ylim=c(1,5)))
with(pdat, lines(mass[i], estimate[i], col=cl))
with(pdat, lines(mass[i], estimate[i]-se[i], lty=2, col=cl))
with(pdat, lines(mass[i], estimate[i]+se[i], lty=2, col=cl))
##################################
#GROUPED DATA
##################################
dat2 <- subset(dat, !is.na(temp) & !is.na(humid))

cuts <- c(0,1,2,3,4,6,10)
dcat <- cut(dat2$radius, cuts)
distbegin <- cuts[dcat]
distend <- cuts[as.numeric(dcat)+1]
dat2$distance <- dat2$radius
dat2$distcat <- paste(distbegin, distend, sep="-")
mod <- ds(dat2, order=0, key="hr", transect="point", formula=~log(mass))
mod1 <- fitdf(distance~1, dat2, transect="point", key="hr", order=0)
mod2 <- fitdf(distcat~1, dat2, transect="point", key="hr", order=0)
plot(mod1$ddf, pdf=T)
plot(mod2$ddf, pdf=T)
mod1$ddf$par
mod1$edd
mod2$ddf$par
mod2$edd

max(dat2$angle)
fov <- 22.5*pi/180
cuts <- fov * c(0,0.5,1,2.5)
acat <- cut(dat2$angle, cuts)
  angbegin <- cuts[acat]
angend <- cuts[as.numeric(acat)+1]
dat2$anglecat <- paste(angbegin, angend, sep="-")
mod1 <- fitdf(angle~1, dat2)
mod2 <- fitdf(anglecat~1, dat2)
mod3 <- fitdf(anglecat~1, dat2, key="hr")
par(mfrow=c(2,1))
plot(mod1$ddf)
plot(mod2$ddf)
plot(mod3$ddf)
mod1$ddf$par
mod1$edd
mod2$ddf$par
mod2$edd
mod3$ddf$par
mod3$edd

mod4 <- fitdf(anglecat~log(mass), dat2)
mod4$ddf$par
mod4$edd
mod1$ddf$criterion
mod2$ddf$criterion
mod3$ddf$criterion
mod4$ddf$criterion
