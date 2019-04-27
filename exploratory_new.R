library(tidyverse)
#library(RColorBrewer)
library(RCurl)
library(data.table)
library(data.table)
library(spatstat)
library(rgdal)
library(maptools)
library(rjson)
library(broom)
library(geojsonio)

setwd("~/Desktop/Academic/Research/SDSS in R/QSO/R")

QSO_lite <- read.csv("../data/QSO_lite.csv", sep = ",")

z_scores <- read.csv("../data/scores_with_z.csv", sep = ",")
z_scores <- z_scores %>% filter(rating != 0)

# omit the code below as QSO is so large
#QSO <- fread("../data/QSO.csv", sep = ",")

##### read in blue stars
mag_range <- 19.0
mag_delta <- 0.2
g_r_difference <- 0.5
targetSqlQuery = paste("select ObjID, ra, dec from star where ra between ",                        0, " and ", 360, " and dec between ",  -10, " and ",
                       70, " and  r between ", mag_range, " and ", 
                       mag_range + mag_delta," and g-r < ", g_r_difference,
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
AND (((flags & 0x100000000000) = 0) or (flags & 0x1000) = 0)", sep = " ")

targetSqlQuery = gsub(pattern = "\n",replacement = " ", x = targetSqlQuery)
urlBase = "http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/SqlSearch?"
blue_stars = getForm(urlBase, cmd = targetSqlQuery, format = "csv")
blue_stars = read.table(text = blue_stars, header=TRUE, sep=",", dec=".", comment.char="#")
##### end read of blue stars #####


##### Making Spatstat Object of Rated QSO's ####

ow <- as.owin(readOGR("./owin/sdss_owin.shp"))
qso_pp <-  with(z_scores, ppp(ra, dec, ow, marks = Z))
qso_pp <- rescale(qso_pp, unitname = "degrees")
qso_diggle <- bw.diggle(qso_pp, hmax = 0.1)
z_scores$density <- density(qso_pp, sigma = qso_diggle, diggle = T, at = "points")
marks(qso_pp) <- data.frame(Z = z_scores$Z, 
                            background = z_scores$density, 
                            rating = z_scores$rating)

##########

q <- geojson_read("../data/mw1.json", what = "sp")
milky_way <- tidy(q)

#plot of 1) high score quasars, 2) density of scored quasars, 3) milky way
z_scores %>% ggplot(aes(x = ra, y = dec)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", bins = 10) + 
  labs(x = "RA", y = "DEC") +
  scale_fill_distiller(palette = 2, direction = 1, guide = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = 'none') +
  geom_point(data = head(z_scores, 100), 
             aes(x = ra, y = dec), col = "firebrick4") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  geom_line(data = milky_way, aes(x = long + 180, y = -lat), 
            col = "salmon", alpha = 0.2)
##### end plot #####

# plot of 1) density of blue stars, 2) high score qso's
blue_stars %>% ggplot(aes(x = ra, y = dec)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", bins = 10) +   
  labs(x = "RA", y = "DEC") +
  scale_fill_gradient(low = "#ffffff", high = "#67AECB", guide = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = 'none') +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  geom_point(aes(x = ra, y = dec), data = head(z_scores, 100))
###### end plot ######

simplify_QSO <- function(QSO, x) {
  QSO <- QSO %>%
    select(SDSS_NAME, RA, DEC, Z, Z_ERR,
           PSFMAG.1, PSFMAG.2, PSFMAG.3, PSFMAG.4, PSFMAG.5,
           PLATE, MJD, FIBERID)
  QSO$ra <- round(QSO$RA, x)
  QSO$dec <- round(QSO$DEC, x)
  QSO
}

QSO_lite <- simplify_QSO(QSO_lite, 5)

get_spectrum <- function(object = QSO_lite) {
  plate <- object$PLATE
  mjd <- object$MJD
  fibre <- object$FIBERID
  url =  paste0("https://dr14.sdss.org/optical/spectrum/view/data/format=csv/spec=full?plateid=",
                plate, "&mjd=", mjd, "&fiberid=",fibre)
  read.csv(file = url) %>% 
    filter(between(Wavelength, 5500, 7000))
  # head(x)
  # x %>% ggplot(aes(x = Wavelength, y = BestFit)) + geom_line()
}


find_star_match <- function(QSO, ra = 10, dec = 00) {
  g <- QSO$PSFMAG.2
  r <- QSO$PSFMAG.3
  i <- QSO$PSFMAG.4
  mySqlQuery <- paste("select ObjID, specObjID, ra, dec, psfmag_u, psfmag_g, psfmag_r, psfmag_i, psfmag_z
  from star
  where
  ra between (", ra-50, ") and (", ra+50, ")
  AND dec between (", dec-50, ") and (", dec+50, ")
  AND specObjID != 0
  AND (psfmag_g - psfmag_r) between (", g-r-0.005, ") and (", g-r+0.005, ")
  AND (psfmag_r - psfmag_i) between (", r-i-0.005, ") and (", r-i+0.005, ")
  -- AND psfmagerr_u < 0.05 and psfmagerr_g < 0.05 and psfmagerr_r < 0.05 and psfmagerr_i < 0.05 and psfmagerr_z < 0.05
  AND type = 6
  AND ((flags & 0x10000000) != 0)       -- detected in BINNED1
  AND ((flags & 0x8100000c00a4) = 0)     -- not EDGE, NOPROFILE, PEAKCENTER,
  -- NOTCHECKED, PSF_FLUX_INTERP, SATURATED, or BAD_COUNTS_ERROR
  AND ((flags & 0x400000000000) = 0) 
  -- not DEBLEND_NOPEAK or small PSF error
  -- (substitute psfmagerr in other band as appropriate)
  AND (((flags & 0x100000000000) = 0) or (flags & 0x1000) = 0)", sep=" ")
  mySqlQuery <- gsub(pattern="\n",replacement=" ",x=mySqlQuery)
  urlBase <- "http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/SqlSearch?"
  X <- getForm(urlBase, cmd = mySqlQuery, format = "csv")
  X <- read.table(text=X, header=TRUE, sep=",", dec=".", comment.char="#")
  if(dim(X)[1]==0) X[1,]=c("1237669703194771471", 
                           "2603176579236390912", 
                           7.781164150, 8.125156825, 
                           16.28, 15.16, 14.77, 14.65, 14.58)
  X$ObjID <- as.character(X$ObjID)
  X$specObjID <- as.character(X$specObjID)
  X
}


find_star_fibers <- function(specObjID) {
  mySqlQuery <- paste("select PLATE, MJD, FIBERID, class, ra, dec
                      from specObj
                      where specObjID = ", specObjID)
  mySqlQuery <- gsub(pattern="\n",replacement=" ",x=mySqlQuery)
  urlBase <- "http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/SqlSearch?"
  X <- getForm(urlBase, cmd = mySqlQuery, format = "csv")
  read.table(text=X, header = TRUE, sep = ",", dec = ".", comment.char = "#")
}

QSO_correlation <- function(index=1, QSO = QSO_lite) {
  z <- find_star_match(QSO[index,])
  if(z$specObjID[1] == "2603176579236390912") {
    return(c(-999, 0, 0, 0, "ERROR", 0, 0))
  } else {
    z1 <- map(z$specObjID, find_star_fibers) %>% 
      rbindlist() %>% 
      filter(class == "STAR") %>% 
      head(1)
    QSO_spectrum <- get_spectrum(object = QSO_lite[index,])
    reference_spectrum <- get_spectrum(z1)
    bind_cols(cor=cor(QSO_spectrum$BestFit, 
                      reference_spectrum$BestFit), z1)
  }
}


z_scores %>% ggplot(aes(x = ra, y = dec)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", bins = 10) +   
  labs(x = "RA", y = "DEC") +
  scale_fill_distiller(palette = 2, direction = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = 'none') +
  geom_point(data = head(z_scores, 100), 
             aes(x = ra, y = dec), col = "firebrick4") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ow <- as.owin(readOGR("../owin/sdss_owin.shp"))
qso_pp <-  ppp(z_scores$ra, z_scores$dec, ow, marks = z_scores$Z)
qso_pp <- rescale(qso_pp, unitname = "degrees")
qso_diggle <- bw.diggle(qso_pp, hmax = 0.1)
z_scores$density <- density(qso_pp, sigma = qso_diggle, diggle = T, at = "points")
marks(qso_pp) <- data.frame(Z = z_scores$Z, background = z_scores$density)

z_scores[1:100,c(3, 4, 23)] %>% 
  filter(between(ra, 100, 150) & between(dec, 20, 60)) %>% 
  arrange(ra)