library(ggplot2)
library(ggrepel)
library(ggpubr)
library(readr)
library(cowplot)
library(dplyr)
library(tidyr)
library(matrixStats)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)
library(stringr)
library(png)
library(ggbreak)
library(corrplot)
library(ggcorrplot)
library(scales)
library(tibble)
library(MASS)
library(ggpmisc)
library(psych)
library(ggforce)
library(ggtext)
library(hexbin)
library(ggh4x)
library(survminer)
library(ggrastr)

ggplot_theme <- function() {

  theme_classic() +
  theme(

    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.title.align = 0.5,
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,0,0,0),
    legend.text = element_text(colour = "black", size = 8),

    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(colour = "black", size = 12, hjust = 0.5),

    plot.title = element_text(size = 14, hjust = 0.5),

    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 12)
  )
} 

ggplot_theme_bw <- function() {

  theme_bw() +
  theme(

    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.title.align = 0.5,
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-5,0,0,0),
    legend.text = element_text(colour = "black", size = 8),

    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(colour = "black", size = 12, hjust = 0.5),

    plot.title = element_text(size = 14, hjust = 0.5),
    strip.background = element_blank(),

    strip.text = element_text(colour = "black", size = 12)
  )
}
