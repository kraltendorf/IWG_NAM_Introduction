# Project: IWG_NAM_Introduction
# Analysis - LD Decay Calculation
# Authors: Kayla R. Altendorf & Garett Heineck
# Date: 01/31/2020

# load required packages
library("vcfR")
library("genetics")
library("dplyr")
library("tidyverse")
library("segmented")
library("ggpointdensity")
library("writexl")
library("readxl")
library("MALDIquant")

# location of github directory
dir <- "/Users/kayla.altendorf/OneDrive - USDA/Publications/IWG NAM Introduction/Scripts for Github/"

# script we're on
script <- c("Linkage Disequillibrium")

# set path to save output
out_path <- paste(dir, script, "/output", sep = "")
folder.dir <- paste(dir, script, "/output/breakpoint_analysis", sep = "/")
dir.create(folder.dir)

# load data from Calculating Pairwise LD script
LD.decay.dat <- read.table(paste(dir, script, "/output/ld_by_family_by_chrom.txt", sep = "/"),
                          header = T)
LD.decay.dat$family <- as.factor(LD.decay.dat$family)

# set ggplot2 theme
LD_theme<- function(base_size = 16) {
  theme_minimal(base_size = base_size) %+replace%
    theme(strip.background = element_rect(fill = "grey85", color = "black", linetype = 1),
          legend.background =  element_rect(fill = "white", linetype = 0),
          legend.position = "bottom",
          strip.text.x = element_text(size = 20),
          strip.background.x = element_rect(fill = "white", colour = NA),
          strip.text.y = element_text(size = 20),
          strip.background.y = element_rect(fill = "white", colour = NA), 
          complete = TRUE)
}

# run through each chromosome, calculating the breakpoints in the graph

######
######
######
#### Chr 1 ####
######
LD.decay.dat_chr1<- filter(LD.decay.dat, 
                           chr == "1") #change here***

ggplot(LD.decay.dat_chr1, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr1<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr1)
mod.seg.chr1<- segmented(mod.glm.chr1, 
                      seg.Z = ~ distance_cm, 
                      psi = list(distance_cm = c(10,
                                                 27)))
  
pred.dat.chr1<- data.frame(distance_cm = round(seq(1,
                                          max(LD.decay.dat_chr1$distance_cm),
                                          length.out= 500),1))
pred.dat.chr1$r2<- round(predict(mod.seg.chr1, newdata = pred.dat.chr1),3)
pred.dat.chr1$chr<- rep(paste("1")) #change here***
  
threshold.LD.chr1<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr1, r2<0.21 & r2>0.19)) %>%
                       group_by(threshold) %>%
                       summarise(location.distance_cm = mean(location.distance_cm),
                                 location.r2 = mean(location.r2))
  
  
mod.output.chr1<- list(model_pred        =pred.dat.chr1,
                  model_break_points=data.frame(mod.seg.chr1$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr1)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr1)),
                  LD.threshold      =threshold.LD.chr1)
  
write_xlsx(mod.output.chr1, paste(folder.dir, 
                               paste("chromosome_", "chr1", ".xlsx", sep = ""), 
                               sep = "/"))
######
######
######





######
######
######
##### Chr 2 ####
######
LD.decay.dat_chr2<- filter(LD.decay.dat, 
                           chr == "2") #change here***

ggplot(LD.decay.dat_chr2, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr2<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr2)
mod.seg.chr2<- segmented(mod.glm.chr2, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    50)))

pred.dat.chr2<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr2$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr2$r2<- round(predict(mod.seg.chr2, newdata = pred.dat.chr2),3)
pred.dat.chr2$chr<- rep(paste("2")) #change here***

threshold.LD.chr2<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr2, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr2<- list(model_pred        =pred.dat.chr2,
                  model_break_points=data.frame(mod.seg.chr2$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr2)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr2)),
                  LD.threshold      =threshold.LD.chr2)

write_xlsx(mod.output.chr2, paste(folder.dir, 
                             paste("chromosome_", "chr2", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
##### Chr 3 ####
######
LD.decay.dat_chr3<- filter(LD.decay.dat, 
                           chr == "3") #change here***

ggplot(LD.decay.dat_chr3, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr3<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr3)
mod.seg.chr3<- segmented(mod.glm.chr3, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    20,
                                                    50)))

pred.dat.chr3<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr3$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr3$r2<- round(predict(mod.seg.chr3, newdata = pred.dat.chr3),3)
pred.dat.chr3$chr<- rep(paste("3")) #change here***

threshold.LD.chr3<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr3, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr3<- list(model_pred        =pred.dat.chr3,
                  model_break_points=data.frame(mod.seg.chr3$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr3)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr3)),
                  LD.threshold      =threshold.LD.chr3)

write_xlsx(mod.output.chr3, paste(folder.dir, 
                             paste("chromosome_", "chr3", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 4 ####
######
LD.decay.dat_chr4<- filter(LD.decay.dat, 
                           chr == "4") #change here***

ggplot(LD.decay.dat_chr4, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr4<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr4)
mod.seg.chr4<- segmented(mod.glm.chr4, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    60)))

pred.dat.chr4<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr4$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr4$r2<- round(predict(mod.seg.chr4, newdata = pred.dat.chr4),3)
pred.dat.chr4$chr<- rep(paste("4")) #change here***

threshold.LD.chr4<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr4, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr4<- list(model_pred        =pred.dat.chr4,
                  model_break_points=data.frame(mod.seg.chr4$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr4)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr4)),
                  LD.threshold      =threshold.LD.chr4)

write_xlsx(mod.output.chr4, paste(folder.dir, 
                             paste("chromosome_", "chr4", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
## Chr 5
######
LD.decay.dat_chr5<- filter(LD.decay.dat, 
                           chr == "5") #change here***

ggplot(LD.decay.dat_chr5, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr5<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr5)
mod.seg.chr5<- segmented(mod.glm.chr5, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    50)))

pred.dat.chr5<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr5$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr5$r2<- round(predict(mod.seg.chr5, newdata = pred.dat.chr5),3)
pred.dat.chr5$chr<- rep(paste("5")) #change here***

threshold.LD.chr5<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr5, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr5<- list(model_pred        =pred.dat.chr5,
                  model_break_points=data.frame(mod.seg.chr5$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr5)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr5)),
                  LD.threshold      =threshold.LD.chr5)

write_xlsx(mod.output.chr5, paste(folder.dir, 
                             paste("chromosome_", "chr5", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 6 ####
######
LD.decay.dat_chr6<- filter(LD.decay.dat, 
                           chr == "6") #change here***

ggplot(LD.decay.dat_chr6, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr6<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr6)
mod.seg.chr6<- segmented(mod.glm.chr6, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    50)))

pred.dat.chr6<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr6$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr6$r2<- round(predict(mod.seg.chr6, newdata = pred.dat.chr6),3)
pred.dat.chr6$chr<- rep(paste("6")) #change here***

threshold.LD.chr6<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr6, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr6<- list(model_pred        =pred.dat.chr6,
                  model_break_points=data.frame(mod.seg.chr6$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr6)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr6)),
                  LD.threshold      =threshold.LD.chr6)

write_xlsx(mod.output.chr6, paste(folder.dir, 
                             paste("chromosome_", "chr6", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 7 ####
######
LD.decay.dat_chr7<- filter(LD.decay.dat, 
                           chr == "7") #change here***

ggplot(LD.decay.dat_chr7, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr7<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr7)
mod.seg.chr7<- segmented(mod.glm.chr7, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    50)))

pred.dat.chr7<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr7$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr7$r2<- round(predict(mod.seg.chr7, newdata = pred.dat.chr7),3)
pred.dat.chr7$chr<- rep(paste("7")) #change here***

threshold.LD.chr7<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr7, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr7<- list(model_pred        =pred.dat.chr7,
                  model_break_points=data.frame(mod.seg.chr7$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr7)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr7)),
                  LD.threshold      =threshold.LD.chr7)

write_xlsx(mod.output.chr7, paste(folder.dir, 
                             paste("chromosome_", "chr7", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 8 ####
######
LD.decay.dat_chr8<- filter(LD.decay.dat, 
                           chr == "8") #change here***

ggplot(LD.decay.dat_chr8, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr8<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr8)
mod.seg.chr8<- segmented(mod.glm.chr8, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    50)))

pred.dat.chr8<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr8$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr8$r2<- round(predict(mod.seg.chr8, newdata = pred.dat.chr8),3)
pred.dat.chr8$chr<- rep(paste("8")) #change here***

threshold.LD.chr8<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr8, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr8<- list(model_pred        =pred.dat.chr8,
                  model_break_points=data.frame(mod.seg.chr8$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr8)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr8)),
                  LD.threshold      =threshold.LD.chr8)

write_xlsx(mod.output.chr8, paste(folder.dir, 
                             paste("chromosome_", "chr8", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 9 ####
######
LD.decay.dat_chr9<- filter(LD.decay.dat, 
                           chr == "9") #change here***

ggplot(LD.decay.dat_chr9, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr9<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr9)
mod.seg.chr9<- segmented(mod.glm.chr9, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    40)))

pred.dat.chr9<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr9$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr9$r2<- round(predict(mod.seg.chr9, newdata = pred.dat.chr9),3)
pred.dat.chr9$chr<- rep(paste("9")) #change here***

threshold.LD.chr9<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr9, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr9<- list(model_pred        =pred.dat.chr9,
                  model_break_points=data.frame(mod.seg.chr9$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr9)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr9)),
                  LD.threshold      =threshold.LD.chr9)

write_xlsx(mod.output.chr9, paste(folder.dir, 
                             paste("chromosome_", "chr9", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 10 ####
######
LD.decay.dat_chr10<- filter(LD.decay.dat, 
                           chr == "10") #change here***

ggplot(LD.decay.dat_chr10, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr10<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr10)
mod.seg.chr10<- segmented(mod.glm.chr10, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27)))

pred.dat.chr10<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr10$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr10$r2<- round(predict(mod.seg.chr10, newdata = pred.dat.chr10),3)
pred.dat.chr10$chr<- rep(paste("10")) #change here***

threshold.LD.chr10<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr10, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr10<- list(model_pred        =pred.dat.chr10,
                  model_break_points=data.frame(mod.seg.chr10$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr10)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr10)),
                  LD.threshold      =threshold.LD.chr10)

write_xlsx(mod.output.chr10, paste(folder.dir, 
                             paste("chromosome_", "chr10", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 11 ####
######
LD.decay.dat_chr11<- filter(LD.decay.dat, 
                           chr == "11") #change here***

ggplot(LD.decay.dat_chr11, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr11<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr11)
mod.seg.chr11<- segmented(mod.glm.chr11, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    50)))

pred.dat.chr11<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr11$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr11$r2<- round(predict(mod.seg.chr11, newdata = pred.dat.chr11),3)
pred.dat.chr11$chr<- rep(paste("11")) #change here***

threshold.LD.chr11<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr11, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr11<- list(model_pred        =pred.dat.chr11,
                  model_break_points=data.frame(mod.seg.chr11$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr11)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr11)),
                  LD.threshold      =threshold.LD.chr11)

write_xlsx(mod.output.chr11, paste(folder.dir, 
                             paste("chromosome_", "chr11", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 12 ####
######
LD.decay.dat_chr12<- filter(LD.decay.dat, 
                           chr == "12") #change here***

ggplot(LD.decay.dat_chr12, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr12<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr12)
mod.seg.chr12<- segmented(mod.glm.chr12, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27)))

pred.dat.chr12<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr12$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr12$r2<- round(predict(mod.seg.chr12, newdata = pred.dat.chr12),3)
pred.dat.chr12$chr<- rep(paste("12")) #change here***

threshold.LD.chr12<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr12, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr12<- list(model_pred        =pred.dat.chr12,
                  model_break_points=data.frame(mod.seg.chr12$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr12)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr12)),
                  LD.threshold      =threshold.LD.chr12)

write_xlsx(mod.output.chr12, paste(folder.dir, 
                             paste("chromosome_", "chr12", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 13 ####
######
LD.decay.dat_chr13<- filter(LD.decay.dat, 
                           chr == "13") #change here***

ggplot(LD.decay.dat_chr13, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr13<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr13)
mod.seg.chr13<- segmented(mod.glm.chr13, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    50)))

pred.dat.chr13<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr13$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr13$r2<- round(predict(mod.seg.chr13, newdata = pred.dat.chr13),3)
pred.dat.chr13$chr<- rep(paste("13")) #change here***

threshold.LD.chr13<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr13, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr13<- list(model_pred        =pred.dat.chr13,
                  model_break_points=data.frame(mod.seg.chr13$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr13)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr13)),
                  LD.threshold      =threshold.LD.chr13)

write_xlsx(mod.output.chr13, paste(folder.dir, 
                             paste("chromosome_", "chr13", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 14 ####
######
LD.decay.dat_chr14<- filter(LD.decay.dat, 
                           chr == "14") #change here***

ggplot(LD.decay.dat_chr14, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr14<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr14)
mod.seg.chr14<- segmented(mod.glm.chr14, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    50)))

pred.dat.chr14<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr14$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr14$r2<- round(predict(mod.seg.chr14, newdata = pred.dat.chr14),3)
pred.dat.chr14$chr<- rep(paste("14")) #change here***

threshold.LD.chr14<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr14, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr14<- list(model_pred        =pred.dat.chr14,
                  model_break_points=data.frame(mod.seg.chr14$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr14)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr14)),
                  LD.threshold      =threshold.LD.chr14)

write_xlsx(mod.output.chr14, paste(folder.dir, 
                             paste("chromosome_", "chr14", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 15 ####
######
LD.decay.dat_chr15<- filter(LD.decay.dat, 
                           chr == "15") #change here***

ggplot(LD.decay.dat_chr15, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr15<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr15)
mod.seg.chr15<- segmented(mod.glm.chr15, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    40)))

pred.dat.chr15<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr15$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr15$r2<- round(predict(mod.seg.chr15, newdata = pred.dat.chr15),3)
pred.dat.chr15$chr<- rep(paste("15")) #change here***

threshold.LD.chr15<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr15, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr15<- list(model_pred        =pred.dat.chr15,
                  model_break_points=data.frame(mod.seg.chr15$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr15)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr15)),
                  LD.threshold      =threshold.LD.chr15)

write_xlsx(mod.output.chr15, paste(folder.dir, 
                             paste("chromosome_", "chr15", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 16 ####
######
LD.decay.dat_chr16<- filter(LD.decay.dat, 
                           chr == "16") #change here***

ggplot(LD.decay.dat_chr16, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr16<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr16)
mod.seg.chr16<- segmented(mod.glm.chr16, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(15)))

pred.dat.chr16<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr16$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr16$r2<- round(predict(mod.seg.chr16, newdata = pred.dat.chr16),3)
pred.dat.chr16$chr<- rep(paste("16")) #change here***

threshold.LD.chr16<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr16, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr16<- list(model_pred        =pred.dat.chr16,
                  model_break_points=data.frame(mod.seg.chr16$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr16)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr16)),
                  LD.threshold      =threshold.LD.chr16)

write_xlsx(mod.output.chr16, paste(folder.dir, 
                             paste("chromosome_", "chr16", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 17 ####
######
LD.decay.dat_chr17<- filter(LD.decay.dat, 
                           chr == "17") #change here***

ggplot(LD.decay.dat_chr17, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr17<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr17)
mod.seg.chr17<- segmented(mod.glm.chr17, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    40)))

pred.dat.chr17<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr17$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr17$r2<- round(predict(mod.seg.chr17, newdata = pred.dat.chr17),3)
pred.dat.chr17$chr<- rep(paste("17")) #change here***

threshold.LD.chr17<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr17, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr17<- list(model_pred        =pred.dat.chr17,
                  model_break_points=data.frame(mod.seg.chr17$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr17)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr17)),
                  LD.threshold      =threshold.LD.chr17)

write_xlsx(mod.output.chr17, paste(folder.dir, 
                             paste("chromosome_", "chr17", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 18 ####
######
LD.decay.dat_chr18<- filter(LD.decay.dat, 
                           chr == "18") #change here***

ggplot(LD.decay.dat_chr18, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr18<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr18)
mod.seg.chr18<- segmented(mod.glm.chr18, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27)))

pred.dat.chr18<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr18$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr18$r2<- round(predict(mod.seg.chr18, newdata = pred.dat.chr18),3)
pred.dat.chr18$chr<- rep(paste("18")) #change here***

threshold.LD.chr18<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr18, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))

mod.output.chr18<- list(model_pred        =pred.dat.chr18,
                  model_break_points=data.frame(mod.seg.chr18$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr18)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr18)),
                  LD.threshold      =threshold.LD.chr18)


write_xlsx(mod.output.chr18, paste(folder.dir, 
                             paste("chromosome_", "chr18", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 19 ####
######
LD.decay.dat_chr19<- filter(LD.decay.dat, 
                           chr == "19") #change here***

ggplot(LD.decay.dat_chr19, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr19<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr19)
mod.seg.chr19<- segmented(mod.glm.chr19, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    60)))

pred.dat.chr19<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr19$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr19$r2<- round(predict(mod.seg.chr19, newdata = pred.dat.chr19),3)
pred.dat.chr19$chr<- rep(paste("19")) #change here***

threshold.LD.chr19<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr19, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr19<- list(model_pred        =pred.dat.chr19,
                  model_break_points=data.frame(mod.seg.chr19$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr19)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr19)),
                  LD.threshold      =threshold.LD.chr19)

write_xlsx(mod.output.chr19, paste(folder.dir, 
                             paste("chromosome_", "chr19", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 20 ####
######
LD.decay.dat_chr20<- filter(LD.decay.dat, 
                           chr == "20") #change here***

ggplot(LD.decay.dat_chr20, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr20<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr20)
mod.seg.chr20<- segmented(mod.glm.chr20, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    40)))

pred.dat.chr20<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr20$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr20$r2<- round(predict(mod.seg.chr20, newdata = pred.dat.chr20),3)
pred.dat.chr20$chr<- rep(paste("20")) #change here***

threshold.LD.chr20<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr20, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr20<- list(model_pred        =pred.dat.chr20,
                  model_break_points=data.frame(mod.seg.chr20$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr20)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr20)),
                  LD.threshold      =threshold.LD.chr20)

write_xlsx(mod.output.chr20, paste(folder.dir, 
                             paste("chromosome_", "chr20", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
#### Chr 21 ####
######
LD.decay.dat_chr21<- filter(LD.decay.dat, 
                           chr == "21") #change here***

ggplot(LD.decay.dat_chr21, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey80", 
                       high="grey20")+
  LD_theme(base_size = 18)

mod.glm.chr21<- glm(r2 ~ distance_cm, data = LD.decay.dat_chr21)
mod.seg.chr21<- segmented(mod.glm.chr21, 
                         seg.Z = ~ distance_cm, 
                         psi = list(distance_cm = c(10,
                                                    27,
                                                    50)))

pred.dat.chr21<- data.frame(distance_cm = round(seq(1,
                                                   max(LD.decay.dat_chr21$distance_cm),
                                                   length.out= 500),1))
pred.dat.chr21$r2<- round(predict(mod.seg.chr21, newdata = pred.dat.chr21),3)
pred.dat.chr21$chr<- rep(paste("21")) #change here***

threshold.LD.chr21<- data.frame(threshold=c(0.2),
                               location=filter(pred.dat.chr21, r2<0.21 & r2>0.19)) %>%
  group_by(threshold) %>%
  summarise(location.distance_cm = mean(location.distance_cm),
            location.r2 = mean(location.r2))


mod.output.chr21<- list(model_pred        =pred.dat.chr21,
                  model_break_points=data.frame(mod.seg.chr21$psi),
                  model_intercepts  =data.frame(intercept(mod.seg.chr21)[[1]]),
                  model_slopes      =data.frame(slope(mod.seg.chr21)),
                  LD.threshold      =threshold.LD.chr21)

write_xlsx(mod.output.chr21, paste(folder.dir, 
                             paste("chromosome_", "chr21", ".xlsx", sep = ""), 
                             sep = "/"))
######
######
######





######
######
######
## Loading LD decay model output
######
LD.files<- list.files(path=folder.dir, full.names = T)
LD.names<- list.files(path=folder.dir, full.names = F)

LD.pred.dat<- data.frame()
for (i in 1:length(LD.files)) {
  temp<- read_excel(paste(LD.files[i]),
             sheet = "model_pred")
  LD.pred.dat<- rbind(temp,LD.pred.dat)
}

LD.threshold<- data.frame()
for (i in 1:length(LD.files)) {
  dat<- read_excel(paste(LD.files[i]),
                    sheet = "LD.threshold")
  temp<- data.frame(distance_cm = dat$location.distance_cm,
                    chr         = sub(".*_chr *(.*?) *.xlsx*", "\\1", LD.names[i]))
  LD.threshold<- rbind(temp,LD.threshold)
}

######
######
######




######
######
######
## Plotting LD decay
######
panel.order<- c("1","4","7","10","13","16","19","2","5","8","11",'14','17',"20","3","6","9","12","15",'18','21')

ggplot(LD.decay.dat, aes(x=distance_cm, y=r2)) +
  geom_pointdensity() +
  scale_color_gradient(low="grey70", 
                       high="grey20")+
  facet_wrap(~factor(chr, levels = panel.order), 
                     scales = "free_x", 
             nrow = 3, 
             ncol = 7) +
  geom_segment(aes(x=0,
                   xend=210,
                   y=.2,
                   yend=.2),
               color="#770026" ,
               size=1,
               linetype="dotted")+
  geom_line(LD.pred.dat, 
            mapping = aes(x=distance_cm, y=r2),
            color="#0072B2",
            size=1.5)+
  geom_text(data = LD.threshold, aes(x = 185, y = .8, 
                                     label = round(distance_cm,1),
                                     size = 5,
                                     color = NULL,group= NULL))+
  labs(x="Pairwise Distance (cM)",
       y=bquote('LD'~(R^2)),
       color="Neighbor Density")+
  LD_theme(base_size = 22)+
  theme(legend.text = element_text(size = 10,
                                   angle= 45),
        legend.title = element_text(size = 18),
        legend.position = "None",
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 15))


# summarise results
median(LD.threshold$distance_cm)
hist(LD.threshold$distance_cm)

# print the figure
ggsave("ld_decay.tiff", 
       plot = last_plot(), 
       path = paste(dir, script, "/output", sep = ""),
       device = "tiff",
       scale = 1, 
       width = 15, 
       height = 10, 
       units = "in")

######
######
######
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
######
######
######


