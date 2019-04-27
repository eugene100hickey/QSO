setwd("~/Desktop/Academic/Research/SDSS in R/QSO")

library(ggplot2)
library(dplyr)
library(RCurl)
library(patchwork)
library(spatstat)

#read and clean up data
df <- read.csv("data/quasar_catalogue.csv", sep=",")
names(df)[4:8] = c("u_band", "g_band", 
                   "r_band", "i_band", "z_band" )

df <- df %>% filter(r_band > 10)
df <- df %>% filter(rating > 0.2)

df %>% ggplot(aes(x=rating)) +
  geom_histogram(col="navyblue", bins=200) + 
  xlim(0,8) +
  theme_minimal()

#plot of rating vs r-band magnitude with loess line fit
df %>% ggplot(aes(x=r_band, y=rating)) + 
  geom_point(shape=1, alpha=0.2, colour="cyan") + 
  geom_smooth(method="loess") +
  theme_minimal()

#sky distribution of quasars, top ranked quasers in red
g1 <- ggplot(data=df, aes(x=ra, y=dec, fill=-rating, col=-rating)) +
  geom_point(shape=21, alpha=0.2, size=0.1, show.legend=F) +
  geom_point(data=df[1:100,], aes(x=ra, y=dec), 
             col="darkred", size=0.2, show.legend=F) +
  theme_minimal()
g2 <- df[1:200, ] %>% ggplot(aes(x=ra, y=dec, fill=-rating, col=-rating)) +
  geom_point(shape=21, alpha=0.8, size=0.5) +
  theme_minimal()

g1

#g-r vs r-i density plot
df %>% mutate(g_r = g_band-r_band, r_i = r_band-i_band) %>% 
  ggplot(aes(x=r_i, y=g_r)) +
  stat_density_2d(aes(fill = ..level..),
                  geom = "polygon", bins=15) +   
  labs(x="r-i magnitude difference", 
       y="g-r magnitude difference") +
  theme_minimal()




#some code to try and download ancilliary data from SDSS
epsilon <- 0.01
ra <- 313.0512
dec <- -0.445892

targetSqlQuery = paste("select ObjID, ra, dec from star where ra between ", 
                       0, " and ", 
                       360, " and dec between ", 
                       -10, " and ", 70, " and  r between ", 19.4, " and ", 19.6,
                       " and g-r < 0.5", 
"--and psfmagerr_u < 0.05 and 
psfmagerr_g < 0.05 and 
psfmagerr_r < 0.05 and 
psfmagerr_i < 0.05 and 
psfmagerr_z < 0.05
and ((flags & 0x10000000) != 0)       -- detected in BINNED1
AND ((flags & 0x8100000c00a4) = 0)     -- not EDGE, NOPROFILE, PEAKCENTER,
-- NOTCHECKED, PSF_FLUX_INTERP, SATURATED, or BAD_COUNTS_ERROR
AND ((flags & 0x400000000000) = 0) 
-- not DEBLEND_NOPEAK or small PSF error
-- (substitute psfmagerr in other band as appropriate)
AND (((flags & 0x100000000000) = 0) or (flags & 0x1000) = 0)", sep=" ")

targetSqlQuery = gsub(pattern="\n",replacement=" ",x=targetSqlQuery)
urlBase = "http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/SqlSearch?"
X = getForm(urlBase, cmd = targetSqlQuery, format = "csv")
X = read.table(text=X, header=TRUE, sep=",", dec=".", comment.char="#")
