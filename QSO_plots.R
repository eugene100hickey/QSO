library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(RColorBrewer)
library(RCurl)

setwd("~/Desktop/Academic/Research/SDSS in R/QSO/R")

QSO_lite <- read.csv("../data/QSO_lite.csv", sep = ",")

# scores <- read.csv("../data/quasar_catalogue.csv", sep = ",")
# scores$ra <- round(scores$ra, 4)
# scores$dec <- round(scores$dec, 4)

z_scores <- read.csv("../data/scores_with_z.csv", sep = ",")


commonTheme <- list(
  labs(
    color = "Z",
    fill = "Z",
    x = "RA",
    y = "DEC"
  ),
  theme_bw(),
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1)
  )
)

#histogram of Z values
QSO_lite %>% ggplot(aes(x = Z)) +
  geom_histogram(bins = 100, fill = "firebrick4") +
  xlim(0, 4.5) +
  ggtitle("Frequency of Z scores from all QSO's of SDSS")

#plot of colour differences
QSO_lite %>%
  mutate(deltaR = PSFMAG.2 - PSFMAG.3,
         deltaI = PSFMAG.3 - PSFMAG.4) %>%
  ggplot(aes(x = deltaI, y = deltaR)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", bins = 15) +
  labs(x = "r-i magnitude difference",
       y = "g-r magnitude difference") +
  theme_minimal()

#plot of colour differences at different Z's
QSO_lite %>%
  mutate(Z_class = cut_interval(QSO_lite$Z, length = 0.5)) %>%
  filter(Z < 4) %>%
  mutate(deltaR = PSFMAG.2 - PSFMAG.3,
         deltaI = PSFMAG.3 - PSFMAG.4) %>%
  ggplot(aes(x = deltaI, y = deltaR)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", bins = 15) +
  labs(x = "r-i magnitude difference",
       y = "g-r magnitude difference") +
  facet_wrap( ~ Z_class) +
  xlim(-0.25, 0.75) +
  ylim(-0.25, 2) +
  theme_minimal()

#RA, DEC plot
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

QSO_lite %>%
  filter(Z > 3) %>%
  ggplot(aes(
    x = RA,
    y = DEC,
    col = Z,
    fill = Z
  )) +
  geom_point(alpha = 0.2, aes(size = Z)) +
  scale_colour_gradientn(colours = myPalette(100),
                         values = (1:100) / 100) +
  geom_density2d(aes(colour = ..level..)) +
  commonTheme


QSO_colour_diff <- QSO_lite %>%
  select(Z, PSFMAG.1, PSFMAG.2, PSFMAG.3, PSFMAG.4, PSFMAG.5) %>%
  na.omit %>%
  filter(Z < 4) %>%
  mutate(deltaG = PSFMAG.1 - PSFMAG.2,
         deltaR = PSFMAG.2 - PSFMAG.3,
         deltaI = PSFMAG.3 - PSFMAG.4,
         deltaZ = PSFMAG.4 - PSFMAG.5) %>%
  select(Z, deltaG, deltaR, deltaI, deltaZ)

wss <- function(k) {
  kmeans(QSO_colour_diff, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares")

km_clust <- QSO_colour_diff %>%
  select(-Z) %>%
  kmeans(centers = 4, nstart = 25)


cluster <- QSO_colour_diff %>%
  cbind(cluster = km_clust$cluster)

cluster %>%
  ggplot(aes(x = Z, color = as.factor(cluster),
             fill = as.factor(cluster))) +
  geom_density(position = "identity", alpha = 0.6) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.6)


cluster %>%
  group_by(cluster) %>%
  summarise(mean_deltaG = mean(deltaG),
            mean_deltaR = mean(deltaR),
            mean_deltaI = mean(deltaI),
            mean_deltaZ = mean(deltaZ))


simplyfy_QSO <- function(QSO, x) {
  QSO <- QSO %>%
    select(SDSS_NAME, RA, DEC, Z, Z_ERR,
           PSFMAG.1, PSFMAG.2, PSFMAG.3, PSFMAG.4, PSFMAG.5,
           PLATE, MJD, FIBERID)
  QSO$ra <- round(QSO$RA, x)
  QSO$dec <- round(QSO$DEC, x)
  QSO
}

QSO_lite <- simplyfy_QSO(QSO_lite, 5)



QSO_lite %>% ggplot(aes(x=ra, y=dec)) + 
  geom_point(alpha=0.1, col = "lightgreen", size=0.1) + 
  geom_point(aes(x=ra, y=dec), data=head(QSO_lite, 100))




