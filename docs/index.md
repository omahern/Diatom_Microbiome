---
title: "CoGrowth"
output: 
  html_document:
    keep_md: true
    toc: true
    toc_float: true
    toc_depth: 6
    code_folding: hide
    number_sections: false
    theme: lumen

knit: (function(input_file, encoding) {
  out_dir <- 'docs';
  rmarkdown::render(input_file,
 encoding=encoding,
 output_file=file.path(dirname(input_file), out_dir, 'index.html'))})
---
# workplace setup

``` r
library(Hmisc)
library(vegan)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(phyloseq)
library(ape)
library(ggpubr)
library(microbiome)
library(ape)
library(rstatix)
library(ggpubr)
library(MicEco)
```


# Physiology
## Growth Rates
### Bacteria 

``` r
read<-read.csv(file="datafiles/Bact_Growth_rates.csv",
               header=T,row.names=1)

boxplot(read$Growth~read$Geno)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/bact-1-1.png)<!-- -->

``` r
a=aov(read$Growth~read$Geno)
summary(a)
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)  
## read$Geno    6 0.4072 0.06787   3.938 0.0277 *
## Residuals   10 0.1723 0.01723                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 4 observations deleted due to missingness
```

``` r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/bact-1-2.png)<!-- -->

``` r
boxplot(read$Growth~read$Pop)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/bact-1-3.png)<!-- -->

``` r
a=aov(read$Growth~read$Pop)
summary(a)
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)  
## read$Pop     2 0.1994 0.09969   3.671 0.0523 .
## Residuals   14 0.3802 0.02716                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 4 observations deleted due to missingness
```

``` r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/bact-1-4.png)<!-- -->

### Diatoms 


### used 11 February 


``` r
read<-read.csv(file="datafiles/TR_Growth_rates_Feb16.csv",
               header=T,row.names=1)

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", " "))

subset=subset(read, treatment!="RD4_YE5" & treatment!="RA5_RD4" & treatment!="RD4_NB7")

 
ggplot1=ggplot(subset, aes(x=strain, y=growth_rate, 
                   fill=culture, 
                   color=strain)) +
  geom_boxplot(alpha=.4, linewidth=0.7) +
  theme_classic() + 
  scale_fill_manual(values=c('white','gray'))+

    scale_color_manual(values=c('#4D9DE0','#E15554', '#3BB273',"#7768AE")) +
  theme(plot.margin=unit(c(0,0,0,0.25),"cm"))+
# theme(legend.position="bottom",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 
   theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 

  ylab(expression("Growth rate " * mu *  ~day^-1)) +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  scale_x_discrete(name='') +coord_flip() +
  stat_compare_means(aes(x=strain, y=growth_rate,
                 color=factor(culture)),paired=T,  label = "p.signif", method='anova', hide.ns=TRUE)
ggplot1
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fcm-growth-2-1.png)<!-- -->

#### with fv.fm

``` r
read<-read.csv(file="datafiles/TR_Growth_rates_Feb16.csv",
               header=T,row.names=1)
dim(read)
```

```
## [1] 24 12
```

``` r
head(read)
```

```
##    treatment isolate treat_num growth_rate         names     id population
## 13       NB7     NB7         5    1.310704  NB7 18 Xenic  xenic       PopC
## 14       NB7     NB7         5    1.497809  NB7 18 Xenic  xenic       PopC
## 15       NB7     NB7         5    1.470415  NB7 18 Xenic  xenic       PopC
## 16    NB7_AX     NB7         6    1.474217 NB7 18 Axenic axenic       PopC
## 17    NB7_AX     NB7         6    1.484045 NB7 18 Axenic axenic       PopC
## 18    NB7_AX     NB7         6    1.222048 NB7 18 Axenic axenic       PopC
##    isolate2 culture strain strain_26     fv
## 13      NB7   xenic NbQ-B7        C1 0.5092
## 14      NB7   xenic NbQ-B7        C1 0.5031
## 15      NB7   xenic NbQ-B7        C1 0.4929
## 16      NB7  axenic NbQ-B7        C1 0.4952
## 17      NB7  axenic NbQ-B7        C1 0.4772
## 18      NB7  axenic NbQ-B7        C1 0.4678
```

``` r
read2=subset(read, strain!="NB4" & strain !="NB6" & strain !="RA4" & strain !="YC5")

gplot2=ggplot(read, 
             aes(x=strain, y=fv,
                 color=factor(strain), fill=factor(culture)))+
  geom_boxplot(alpha=.4, linewidth=0.7) +
  theme_classic() + 
  scale_fill_manual(values=c('white','gray'))+
 # scale_color_manual(values=c('#003f5c','#7a5195', '#ef5675',"#ffa600")) +
      scale_color_manual(values=c('#4D9DE0','#E15554', '#3BB273',"#7768AE")) +

    theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 

 # scale_y_continuous(name='Growth Rate')+
  ylab("(Fm/Fo)/Fm") +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  scale_x_discrete(name='')+ coord_flip()
cowplot::plot_grid(ggplot1, gplot2, 
                   rel_heights = c(1,1),
                   nrow=2,
                   align='h',
                   axis='l',
                   byrow=TRUE, labels='auto')
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-growth-a-1.png)<!-- -->

``` r
       #           rel_widths = c(1.7,1))
```


#### fv/fm with sign

``` r
read<-read.csv(file="datafiles/TR_Growth_rates_Feb16.csv",
               header=T,row.names=1)
dim(read)
```

```
## [1] 24 12
```

``` r
head(read)
```

```
##    treatment isolate treat_num growth_rate         names     id population
## 13       NB7     NB7         5    1.310704  NB7 18 Xenic  xenic       PopC
## 14       NB7     NB7         5    1.497809  NB7 18 Xenic  xenic       PopC
## 15       NB7     NB7         5    1.470415  NB7 18 Xenic  xenic       PopC
## 16    NB7_AX     NB7         6    1.474217 NB7 18 Axenic axenic       PopC
## 17    NB7_AX     NB7         6    1.484045 NB7 18 Axenic axenic       PopC
## 18    NB7_AX     NB7         6    1.222048 NB7 18 Axenic axenic       PopC
##    isolate2 culture strain strain_26     fv
## 13      NB7   xenic NbQ-B7        C1 0.5092
## 14      NB7   xenic NbQ-B7        C1 0.5031
## 15      NB7   xenic NbQ-B7        C1 0.4929
## 16      NB7  axenic NbQ-B7        C1 0.4952
## 17      NB7  axenic NbQ-B7        C1 0.4772
## 18      NB7  axenic NbQ-B7        C1 0.4678
```

``` r
read2=subset(read, strain!="NB4" & strain !="NB6" & strain !="RA4" & strain !="YC5")

gplot2=ggplot(read, 
             aes(x=strain, y=fv,
                 color=factor(strain), fill=factor(culture)))+
  geom_boxplot(alpha=.4, linewidth=0.7) +
  theme_classic() + 
  scale_fill_manual(values=c('white','gray'))+
 # scale_color_manual(values=c('#003f5c','#7a5195', '#ef5675',"#ffa600")) +
      scale_color_manual(values=c('#4D9DE0','#E15554', '#3BB273',"#7768AE")) +

    theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 

 # scale_y_continuous(name='Growth Rate')+
  ylab("(Fm/Fo)/Fm") +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  scale_x_discrete(name='')+ coord_flip() +
  stat_compare_means(aes(x=strain, y=fv,
                 color=factor(culture)),paired=T,  label = "p.signif", method='anova', hide.ns=TRUE)
    

gplot2
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-growth-b-1.png)<!-- -->

``` r
cowplot::plot_grid(ggplot1, gplot2, 
                   rel_heights = c(1,1),
                   nrow=2,
                   align='h',
                   axis='l',
                   byrow=TRUE, labels='auto')
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-growth-b-2.png)<!-- -->

``` r
       #           rel_widths = c(1.7,1))
```



#### with terminal biomass



``` r
counts=read.csv(file='datafiles/rotula_cellcounts.csv',header=T,row.names=1)

dim(counts)
```

```
## [1] 24  9
```

``` r
head(counts)
```

```
##          id strain bacteria treatment       X0       X2       X3       X4
## 13    NB7-A NbQ-B7    xenic       NB7 424.4707 6599.328 21654.18 35283.54
## 14    NB7-B NbQ-B7    xenic       NB7 344.3819 6983.744 30797.13 37974.73
## 15    NB7-C NbQ-B7    xenic       NB7 304.3375 7876.801 25068.91 40810.26
## 16 NB7-AX_A NbQ-B7   axenic    NB7_AX 344.3819 8752.541 28692.84 61377.68
## 17 NB7-AX_B NbQ-B7   axenic    NB7_AX 240.2664 5518.121 20617.26 50948.51
## 18 NB7-AX_C NbQ-B7   axenic    NB7_AX 344.3819 5013.560 13465.61 40546.63
##          X5
## 13 25264.35
## 14 26216.47
## 15 34929.60
## 16 54533.74
## 17 48081.93
## 18 50911.20
```

``` r
gplot_counts=ggplot(counts, 
             aes(x=strain, y=X5,
                 color=factor(strain), fill=factor(bacteria)))+
  geom_boxplot(alpha=.4, linewidth=0.7) +
  theme_classic() + 
  scale_fill_manual(values=c('white','gray'))+
      scale_color_manual(values=c('#4D9DE0','#E15554', '#3BB273',"#7768AE")) +

    theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 

 # scale_y_continuous(name='Growth Rate')+
  ylab("Terminal Biomass (cells/mL)") +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  scale_x_discrete(name='')+ coord_flip() +
  stat_compare_means(aes(x=strain, y=X5,
                 color=factor(bacteria)),paired=T,  label = "p.signif", method='anova', hide.ns=TRUE)

gplot_counts
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/final-growth-1.png)<!-- -->

``` r
cowplot::plot_grid(ggplot1,  gplot_counts,   gplot2,
                   rel_heights = c(1,1),
                   nrow=1,
                   align='h',
                   axis='l',
                   byrow=TRUE, labels='auto')
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/final-growth-2.png)<!-- -->

``` r
       #           rel_widths = c(1.7,1))
```

##### stats of fv/fm




## Growth Curves

### Diatoms

``` r
read=read.csv(file='datafiles/rotula_cellcounts.csv',header=T,row.names=1)
dim(read)
```

```
## [1] 24  9
```

``` r
treatment=as.factor(read$treatment[1:48])

dim(read)
```

```
## [1] 24  9
```

``` r
growth_data=data.frame(read[1:48,3:7])

# get mean and standard deviation
mean=aggregate(growth_data,by=list(treatment),FUN=mean,na.rm=T)
mean=t(mean[,2:6])
sd=aggregate(growth_data,by=list(treatment),FUN=sd,na.rm=T)
sd=t(sd[,2:6])
# set time
time=c(0,2,3,4,5)

doubles=read[61:69,]
treat_doub=doubles$treatment
doubles_data=doubles[,3:7]

mean_db=aggregate(doubles_data,by=list(treat_doub),FUN=mean,na.rm=T)
mean_db=t(mean_db[,2:6])
sd_db=aggregate(doubles_data,by=list(treat_doub),FUN=sd,na.rm=T)
sd_db=t(sd_db[,2:6])
```


#### Plot Subset of Strains (used)

``` r
read=read.csv(file='datafiles/rotula_cellcounts.csv',header=T)
# dim(read)

growth_data=data.frame(read[,6:10])

# get mean and standard deviation
mean=aggregate(growth_data,by=list(read$treatment),FUN=mean,na.rm=T)
sd=aggregate(growth_data,by=list(read$treatment),FUN=sd,na.rm=T)
# set time
time=c(0,2,3,4,5)
#     scale_color_manual(values=c('#4D9DE0','#E15554', '#3BB273',"#7768AE")) +

{par(mfrow=c(2,2),
     mar=c(5,5,1,1))
  
  mean_ra5=as.numeric(subset(mean, Group.1=="RA5")[,2:6])
  sd_ra5=as.numeric(subset(sd, Group.1=="RA5")[,2:6])
  #RA5
  errbar(x=time,y=mean_ra5,
         yplus=mean_ra5-sd_ra5,
         yminus=mean_ra5+sd_ra5,
         ylim=c(100,50000),log='y',type='b',col='#7768AE', pch=21, bg='gray',errbar.col = '#7768AE',cex=1.3,
         xlab="Time (days)",
         ylab="Diatom Concentration (cells/mL)",xaxt='n',
         yaxt='n')
  mean_ra5=as.numeric(subset(mean, Group.1=="RA5_AX")[,2:6])
  sd_ra5=as.numeric(subset(sd, Group.1=="RA5_AX")[,2:6])
  errbar(x=time,y=mean_ra5,
         yplus=mean_ra5-sd_ra5,
         yminus=mean_ra5+sd_ra5,
         col='#7768AE', pch=21, bg='white',errbar.col = '#7768AE',cex=1.3,
         type='b',add=T)
  text(1,200,"A1",font=2)
  axis(2,at=c(100,1000,12500,50000))
  axis(1, at=c(0:5),label=TRUE)
  legend("bottomright", legend=c("xenic",'axenic'), pch=21, pt.bg=c('grey','white'), bty='n')

  
  ### RD4
mean_ra5=as.numeric(subset(mean, Group.1=="RD4")[,2:6])
  sd_ra5=as.numeric(subset(sd, Group.1=="RD4")[,2:6])
  errbar(x=time,y=mean_ra5,
         yplus=mean_ra5-sd_ra5,
         yminus=mean_ra5+sd_ra5,
         ylim=c(100,50000),log='y',type='b',col='#4D9DE0', pch=21, bg='gray',errbar.col = '#4D9DE0',cex=1.3,
         xlab="Time (days)",
         ylab="Diatom Concentration (cells/mL)",xaxt='n',
         yaxt='n')
  mean_ra5=as.numeric(subset(mean, Group.1=="RD4_AX")[,2:6])
  sd_ra5=as.numeric(subset(sd, Group.1=="RD4_AX")[,2:6])
  errbar(x=time,y=mean_ra5,
         yplus=mean_ra5-sd_ra5,
         yminus=mean_ra5+sd_ra5,
         col='#4D9DE0', pch=21, bg='white',errbar.col = '#4D9DE0',cex=1.3,
         type='b',add=T)
  text(1,200,"A2",font=2)
  axis(2,at=c(100,1000,12500,50000))
  axis(1, at=c(0:5),label=TRUE)
  legend("bottomright", legend=c("xenic",'axenic'), pch=21, pt.bg=c('grey','white'), bty='n')
  
  
  
   ### RD4
mean_ra5=as.numeric(subset(mean, Group.1=="YE5")[,2:6])
  sd_ra5=as.numeric(subset(sd, Group.1=="YE5")[,2:6])
  errbar(x=time,y=mean_ra5,
         yplus=mean_ra5-sd_ra5,
         yminus=mean_ra5+sd_ra5,
         ylim=c(100,50000),log='y',type='b',col='#E15554', pch=21, bg='gray',errbar.col = '#E15554',cex=1.3,
         xlab="Time (days)",
         ylab="Diatom Concentration (cells/mL)",xaxt='n',
         yaxt='n')
  mean_ra5=as.numeric(subset(mean, Group.1=="YE5_AX")[,2:6])
  sd_ra5=as.numeric(subset(sd, Group.1=="YE5_AX")[,2:6])
  errbar(x=time,y=mean_ra5,
         yplus=mean_ra5-sd_ra5,
         yminus=mean_ra5+sd_ra5,
         col='#E15554', pch=21, bg='white',errbar.col = '#E15554',cex=1.3,
         type='b',add=T)
  text(1,200,"B1",font=2)
  axis(2,at=c(100,1000,12500,50000))
  axis(1, at=c(0:5),label=TRUE)
  legend("bottomright", legend=c("xenic",'axenic'), pch=21, pt.bg=c('grey','white'), bty='n')

  
   ### RD4
mean_ra5=as.numeric(subset(mean, Group.1=="NB7")[,2:6])
  sd_ra5=as.numeric(subset(sd, Group.1=="NB7")[,2:6])
  errbar(x=time,y=mean_ra5,
         yplus=mean_ra5-sd_ra5,
         yminus=mean_ra5+sd_ra5,
         ylim=c(100,50000),log='y',type='b',col='#3BB273', pch=21, bg='gray',errbar.col = '#3BB273',cex=1.3,
         xlab="Time (days)",
         ylab="Diatom Concentration (cells/mL)",xaxt='n',
         yaxt='n')
  mean_ra5=as.numeric(subset(mean, Group.1=="NB7_AX")[,2:6])
  sd_ra5=as.numeric(subset(sd, Group.1=="NB7_AX")[,2:6])
  errbar(x=time,y=mean_ra5,
         yplus=mean_ra5-sd_ra5,
         yminus=mean_ra5+sd_ra5,
         col='#3BB273', pch=21, bg='white',errbar.col = '#3BB273',cex=1.3,
         type='b',add=T)
  text(1,200,"C1",font=2)
  axis(2,at=c(100,1000,12500,50000))
  axis(1, at=c(0:5),label=TRUE)
  legend("bottomright", legend=c("xenic",'axenic'), pch=21, pt.bg=c('grey','white'), bty='n')
  
  }
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/curves-5-1.png)<!-- -->


``` r
read=read.csv(file="datafiles/bact_growth.csv",header=T)
bact1=subset(read,treat !="Bact" )


bact11=mean(bact1$X3)
bact12=mean(bact1$X4)

read=read.csv(file='datafiles/rotula_cellcounts.csv',header=T)
xen=subset(read, bacteria=="xenic")

xen
```

```
##     X    id  strain bacteria treatment       X0       X2        X3       X4
## 1  13 NB7-A  NbQ-B7    xenic       NB7 424.4707 6599.328 21654.183 35283.54
## 2  14 NB7-B  NbQ-B7    xenic       NB7 344.3819 6983.744 30797.130 37974.73
## 3  15 NB7-C  NbQ-B7    xenic       NB7 304.3375 7876.801 25068.906 40810.26
## 7  25 RA5-A  NbO-A5    xenic       RA5 376.4174 2763.064  7105.930 12213.56
## 8  26 RA5-B  NbO-A5    xenic       RA5 344.3820 2610.896  5021.308 26154.83
## 9  27 RA5-C  NbO-A5    xenic       RA5 400.4441 2490.762  6737.420 12061.38
## 13 31 RD4-A  NbO-D4    xenic       RD4 352.3908 2923.242  8612.966 16314.12
## 14 32 RD4-B  NbO-D4    xenic       RD4 424.4707 3043.375  7253.557 19598.40
## 15 33 RD4-C  NbO-D4    xenic       RD4 512.5684 3083.419  9327.747 24493.88
## 19 43 YE5-A NbP-YE5    xenic       YE5 400.4441 3099.437 11363.279 26856.16
## 20 44 YE5-B NbP-YE5    xenic       YE5 344.3819 3195.544 10312.712 25461.12
## 21 45 YE5-C NbP-YE5    xenic       YE5 536.5950 4252.716 13724.118 33113.61
##          X5
## 1  25264.35
## 2  26216.47
## 3  34929.60
## 7  27326.52
## 8  23303.23
## 9  25340.11
## 13 27834.76
## 14       NA
## 15 30441.76
## 19 32990.39
## 20 34573.20
## 21 42021.67
```

``` r
mean1=mean(xen$X3)
mean2=mean(xen$X4)

bact11/mean1
```

```
## [1] 85.38914
```

``` r
bact12/mean2
```

```
## [1] 28.91152
```


## FCM Spectra

#### with lables 


``` r
## add diatom growht rates first 
library(ggridges)
read=read.csv(file='datafiles/fcm_rotula.csv',header=T)
data1=subset(read, Temp !="22")

data1$names=paste(data1$Diatom, data1$Temp,data1$Bact)
data1$naems2=paste(data1$strain_2)



plot3=ggplot(data1, aes(RED,y=naems2, fill=factor(Treat))) +
  geom_density_ridges(alpha=0.5) +
  theme_classic(base_size=12)+
  scale_fill_manual(values=c('#631879','white','#3B4992','white', '#EE0000','white',"#008B45",'white')) +
  theme(plot.margin=unit(c(0,0,0,0),"cm"))+
  # facet_grid(~Pop) +
  theme(legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_continuous(name='692/40 mm', limits=quantile(data1$RED, c(0.01, 0.9999))) + #+facet_wrap(~Pop, scales="free")
  scale_y_discrete(name='')
plot3
```

```
## Picking joint bandwidth of 0.017
```

```
## Warning: Removed 367 rows containing non-finite outside the scale range
## (`stat_density_ridges()`).
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

``` r
ggplot(data1, aes(RED,y=naems2, fill=factor(Treat))) + 
  geom_violin(alpha=0.7, trim=TRUE, scale='area', draw_quantiles  = c(0.25, 0.5, 0.75)) +
  theme_classic(base_size=12)+
   scale_fill_manual(values=c('#631879','white','#3B4992','white', '#EE0000','white',"#008B45",'white')) +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) +
  scale_x_continuous(name='692/40 mm', limits=quantile(data1$RED, c(0.01, 0.9999))) +   coord_flip() 
```

```
## Warning: Removed 367 rows containing non-finite outside the scale range
## (`stat_ydensity()`).
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

``` r
ggplot(data1, aes(RED,y=naems2, fill=factor(Treat))) + 
  geom_violin(alpha=0.7, trim=TRUE, draw_quantiles  = c( 0.5)) +
  theme_classic(base_size=12)+
   scale_fill_manual(values=c('#631879','white','#3B4992','white', '#EE0000','white',"#008B45",'white')) +
  scale_x_continuous(name='692/40 mm', limits=quantile(data1$RED, c(0.01, 0.9999))) +   coord_flip() 
```

```
## Warning: Removed 367 rows containing non-finite outside the scale range
## (`stat_ydensity()`).
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-2-3.png)<!-- -->

``` r
cowplot::plot_grid(plot3,ncol=3, labels='auto')
```

```
## Picking joint bandwidth of 0.017
```

```
## Warning: Removed 367 rows containing non-finite outside the scale range
## (`stat_density_ridges()`).
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-2-4.png)<!-- -->

``` r
levels(as.factor(data1$Treat))
```

```
## [1] "NbO-A5"       "NbO-A5-AX"    "NbO-D4_18"    "NbO-D4_18_AX" "NbP-YE5"     
## [6] "NbP-YE5-AX"   "NbQ-B7-18"    "NbQ-B7-18-AX"
```



#### RED with FIRE 


``` r
## add diatom growht rates first 
library(ggridges)
read=read.csv(file='datafiles/fcm_rotula.csv',header=T)
data1=subset(read, Temp !="22")

data1$names=paste(data1$Diatom, data1$Temp,data1$Bact)
data1$naems2=paste(data1$strain_2)



plot3=ggplot(data1, aes(RED,y=naems2, fill=factor(Treat))) +
  geom_density_ridges(alpha=0.5) +
  theme_classic(base_size=12)+
  scale_fill_manual(values=c('#631879','white','#3B4992','white', '#EE0000','white',"#008B45",'white')) +
  theme(plot.margin=unit(c(0,0,0,0),"cm"))+
  # facet_grid(~Pop) +
  theme(legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_continuous(name='692/40 mm', limits=quantile(data1$RED, c(0.01, 0.9999))) + #+facet_wrap(~Pop, scales="free")
  scale_y_discrete(name='')
plot3
```

```
## Picking joint bandwidth of 0.017
```

```
## Warning: Removed 367 rows containing non-finite outside the scale range
## (`stat_density_ridges()`).
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

``` r
ggplot(data1, aes(RED,y=naems2, fill=factor(Treat))) + 
  geom_violin(alpha=0.7, trim=TRUE, scale='area', draw_quantiles  = c(0.25, 0.5, 0.75)) +
  theme_classic(base_size=12)+
   scale_fill_manual(values=c('#631879','white','#3B4992','white', '#EE0000','white',"#008B45",'white')) +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) +
  scale_x_continuous(name='692/40 mm', limits=quantile(data1$RED, c(0.01, 0.9999))) +   coord_flip() 
```

```
## Warning: Removed 367 rows containing non-finite outside the scale range
## (`stat_ydensity()`).
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

``` r
ggplot(data1, aes(RED,y=naems2, fill=factor(Treat))) + 
  geom_violin(alpha=0.7, trim=TRUE, draw_quantiles  = c( 0.5)) +
  theme_classic(base_size=12)+
   scale_fill_manual(values=c('#631879','white','#3B4992','white', '#EE0000','white',"#008B45",'white')) +
  scale_x_continuous(name='692/40 mm', limits=quantile(data1$RED, c(0.01, 0.9999))) 
```

```
## Warning: Removed 367 rows containing non-finite outside the scale range
## (`stat_ydensity()`).
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-3-3.png)<!-- -->

``` r
cowplot::plot_grid(plot3,ncol=3, labels='auto')
```

```
## Picking joint bandwidth of 0.017
```

```
## Warning: Removed 367 rows containing non-finite outside the scale range
## (`stat_density_ridges()`).
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-3-4.png)<!-- -->

``` r
levels(as.factor(data1$Treat))
```

```
## [1] "NbO-A5"       "NbO-A5-AX"    "NbO-D4_18"    "NbO-D4_18_AX" "NbP-YE5"     
## [6] "NbP-YE5-AX"   "NbQ-B7-18"    "NbQ-B7-18-AX"
```





``` r
read<-read.csv(file="datafiles/TR_Growth_rates_Feb16.csv",
               header=T,row.names=1)
# values_colb=c("#631879", "#3B4992", "#EE0000", "#008B45",
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", " "))
 subset=read
ggplot1=ggplot(subset, aes(x=strain_26, y=growth_rate, 
                   fill=culture, 
                   color=strain_26)) +
  geom_boxplot(alpha=.7) +
  theme_classic(base_size = 14) + 
  scale_fill_manual(values=c('white','gray'))+

    scale_color_manual(values=c('#631879','#3B4992', '#EE0000',"#008B45")) +
  theme(plot.margin=unit(c(0,0,0,0.25),"cm"))+
# theme(legend.position="bottom",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 
   theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1)) + theme(text = element_text(size = 14)) +

  ylab(expression("Growth rate " * mu *  ~day^-1)) +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  scale_x_discrete(name='') +coord_flip() +
  stat_compare_means(aes(x=strain_26, y=growth_rate,
                 color=factor(culture)),paired=T,  label = "p.signif", method='anova', hide.ns=TRUE)
ggplot1
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-growth-c-1.png)<!-- -->

``` r
read<-read.csv(file="datafiles/TR_Growth_rates_Feb16.csv",
               header=T,row.names=1)
dim(read)
```

```
## [1] 24 12
```

``` r
head(read)
```

```
##    treatment isolate treat_num growth_rate         names     id population
## 13       NB7     NB7         5    1.310704  NB7 18 Xenic  xenic       PopC
## 14       NB7     NB7         5    1.497809  NB7 18 Xenic  xenic       PopC
## 15       NB7     NB7         5    1.470415  NB7 18 Xenic  xenic       PopC
## 16    NB7_AX     NB7         6    1.474217 NB7 18 Axenic axenic       PopC
## 17    NB7_AX     NB7         6    1.484045 NB7 18 Axenic axenic       PopC
## 18    NB7_AX     NB7         6    1.222048 NB7 18 Axenic axenic       PopC
##    isolate2 culture strain strain_26     fv
## 13      NB7   xenic NbQ-B7        C1 0.5092
## 14      NB7   xenic NbQ-B7        C1 0.5031
## 15      NB7   xenic NbQ-B7        C1 0.4929
## 16      NB7  axenic NbQ-B7        C1 0.4952
## 17      NB7  axenic NbQ-B7        C1 0.4772
## 18      NB7  axenic NbQ-B7        C1 0.4678
```

``` r
read2=subset(read, strain!="NB4" & strain !="NB6" & strain !="RA4" & strain !="YC5")



gplot2=ggplot(read, 
             aes(x=strain_26, y=fv,
                 color=factor(strain_26), fill=factor(culture)))+
  geom_boxplot(alpha=.7) +
  theme_classic(base_size=14) + 
  scale_fill_manual(values=c('white','gray'))+
 # scale_color_manual(values=c('#003f5c','#7a5195', '#ef5675',"#ffa600")) +
    scale_color_manual(values=c('#631879','#3B4992', '#EE0000',"#008B45")) +
    theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 
 # scale_y_continuous(name='Growth Rate')+
  ylab("Fv/Fm") +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  scale_x_discrete(name='')+ coord_flip() +
  stat_compare_means(aes(x=strain_26, y=fv,
                 color=factor(culture)),paired=T,  label = "p.signif", method='anova', hide.ns=TRUE)
    

gplot2
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-growth-c-2.png)<!-- -->

``` r
gplot22=ggplot(read, 
             aes(x=strain_26, y=fv,
                 color=factor(strain_26), fill=factor(culture)))+
  geom_boxplot(alpha=.7) +
  theme_classic(base_size=14) + 
  scale_fill_manual(values=c('white','gray'))+
 # scale_color_manual(values=c('#003f5c','#7a5195', '#ef5675',"#ffa600")) +
    scale_color_manual(values=c('#631879','#3B4992', '#EE0000',"#008B45")) +
    theme(legend.position="bottom",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 
 # scale_y_continuous(name='Growth Rate')+
  ylab("Fv/Fm") +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  scale_x_discrete(name='')+ coord_flip() +
  stat_compare_means(aes(x=strain_26, y=fv,
                 color=factor(culture)),paired=T,  label = "p.signif", method='anova', hide.ns=TRUE)
  
gplot22
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-growth-c-3.png)<!-- -->

``` r
#ggsave("culture_legend.svg", plot = gplot22, width = 8, height = 5)

gplot3=ggplot(data1, aes(RED,y=naems2, color=factor(strain_2), fill=factor(Bact))) + 
  geom_violin(alpha=0.7, trim=TRUE,draw_quantiles = c(0.5)) +
  theme_classic(base_size=14)+
  ylab("") +
  scale_fill_manual(values=c('white','gray'))+
 # scale_color_manual(values=c('#003f5c','#7a5195', '#ef5675',"#ffa600")) +
    scale_color_manual(values=c('#631879','#3B4992', '#EE0000',"#008B45")) +
  scale_x_continuous(name='692/40 mm', limits=quantile(data1$RED, c(0.01, 0.9999)))  + theme(legend.position = 'none')
gplot3
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-growth-c-4.png)<!-- -->

``` r
data2=data1
trimmed_data <- data2 %>%
  group_by(names) %>%
  filter(
    RED >= quantile(RED, 0.025) &
    RED <= quantile(RED, 0.975)
  ) %>%
  ungroup()


gplot3=ggplot(trimmed_data, aes(RED,y=naems2, color=factor(strain_2), fill=factor(Bact))) + 
  geom_boxplot(alpha=0.7, outliers=FALSE) +
  theme_classic(base_size=14)+
  ylab("") +
  scale_fill_manual(values=c('white','gray'))+
 # scale_color_manual(values=c('#003f5c','#7a5195', '#ef5675',"#ffa600")) +
    scale_color_manual(values=c('#631879','#3B4992', '#EE0000',"#008B45")) +
  scale_x_continuous(name='692/40 mm', limits=quantile(data1$RED, c(0.01, 0.9999)))  + theme(legend.position = 'none')
gplot3
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-growth-c-5.png)<!-- -->

``` r
gplot3=ggplot(trimmed_data, aes(RED,y=naems2, color=factor(strain_2), fill=factor(Bact))) + 
  geom_violin(alpha=0.7, trim=TRUE, draw_quantiles = c(0.5)) +
  theme_classic(base_size=14)+
  ylab("") +
  xlab("692/40 mm") + 
  scale_fill_manual(values=c('white','gray'))+
 # scale_color_manual(values=c('#003f5c','#7a5195', '#ef5675',"#ffa600")) +
    scale_color_manual(values=c('#631879','#3B4992', '#EE0000',"#008B45")) +
   theme(legend.position = 'none')
gplot3
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-growth-c-6.png)<!-- -->

``` r
clams=cowplot::plot_grid(ggplot1, gplot2, gplot3,
                   rel_heights = c(1,1),
                   nrow=1,
                   align='h',
                   axis='l',
                   byrow=TRUE, labels='auto')
#ggsave("filename.svg", plot = clams, width = 8, height = 4)
```


##### stats

``` r
read=read.csv(file='datafiles/fcm_rotula.csv',header=T)

#NbQ-B7
nbq=subset(read, strain=="NbQ-B7")

ks.test(subset(nbq, Bact=="Axenic")$FSC,subset(nbq, Bact=="Xenic")$FSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$FSC and subset(nbq, Bact == "Xenic")$FSC
## D^+ = 0.19175, p-value = 0.0004998
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$SSC,subset(nbq, Bact=="Xenic")$SSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$SSC and subset(nbq, Bact == "Xenic")$SSC
## D^+ = 0.24589, p-value = 0.0004998
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$RED,subset(nbq, Bact=="Xenic")$RED, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$RED and subset(nbq, Bact == "Xenic")$RED
## D^+ = 0.47943, p-value = 0.0004998
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
#NbP-YE5
nbq=subset(read, strain=="NbP-YE5")

ks.test(subset(nbq, Bact=="Axenic")$FSC,subset(nbq, Bact=="Xenic")$FSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$FSC and subset(nbq, Bact == "Xenic")$FSC
## D^+ = 0.0010239, p-value = 0.9945
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$SSC,subset(nbq, Bact=="Xenic")$SSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$SSC and subset(nbq, Bact == "Xenic")$SSC
## D^+ = 0.0030962, p-value = 0.9675
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$RED,subset(nbq, Bact=="Xenic")$RED, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$RED and subset(nbq, Bact == "Xenic")$RED
## D^+ = 0.49622, p-value = 0.0004998
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
#NbO-D4
nbq=subset(read, strain=="NbO-D4")

ks.test(subset(nbq, Bact=="Axenic")$FSC,subset(nbq, Bact=="Xenic")$FSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$FSC and subset(nbq, Bact == "Xenic")$FSC
## D^+ = 0.069753, p-value = 0.0004998
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$SSC,subset(nbq, Bact=="Xenic")$SSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$SSC and subset(nbq, Bact == "Xenic")$SSC
## D^+ = 0.1177, p-value = 0.0004998
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$RED,subset(nbq, Bact=="Xenic")$RED, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$RED and subset(nbq, Bact == "Xenic")$RED
## D^+ = 0.23379, p-value = 0.0004998
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
#NbO-A5
nbq=subset(read, strain=="NbO-A5")

ks.test(subset(nbq, Bact=="Axenic")$FSC,subset(nbq, Bact=="Xenic")$FSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$FSC and subset(nbq, Bact == "Xenic")$FSC
## D^+ = 0.0020021, p-value = 0.9775
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$SSC,subset(nbq, Bact=="Xenic")$SSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$SSC and subset(nbq, Bact == "Xenic")$SSC
## D^+ = 0.0036998, p-value = 0.941
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$RED,subset(nbq, Bact=="Xenic")$RED, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$RED and subset(nbq, Bact == "Xenic")$RED
## D^+ = 0.0048131, p-value = 0.9145
## alternative hypothesis: the CDF of x lies above that of y
```
   
   
##### stats

``` r
read=read.csv(file='datafiles/fcm_rotula.csv',header=T)

#NbQ-B7
nbq=subset(read, strain=="NbQ-B7")

ks.test(subset(nbq, Bact=="Axenic")$FSC,subset(nbq, Bact=="Xenic")$FSC, alternative = 'two.sided', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$FSC and subset(nbq, Bact == "Xenic")$FSC
## D = 0.19175, p-value = 0.0004998
## alternative hypothesis: two-sided
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$SSC,subset(nbq, Bact=="Xenic")$SSC, alternative = 'two.sided', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$SSC and subset(nbq, Bact == "Xenic")$SSC
## D = 0.24589, p-value = 0.0004998
## alternative hypothesis: two-sided
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$RED,subset(nbq, Bact=="Xenic")$RED, alternative = 'two.sided', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$RED and subset(nbq, Bact == "Xenic")$RED
## D = 0.47943, p-value = 0.0004998
## alternative hypothesis: two-sided
```

``` r
#NbP-YE5
nbq=subset(read, strain=="NbP-YE5")

ks.test(subset(nbq, Bact=="Axenic")$FSC,subset(nbq, Bact=="Xenic")$FSC, alternative = 'two.sided', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$FSC and subset(nbq, Bact == "Xenic")$FSC
## D = 0.082073, p-value = 0.0004998
## alternative hypothesis: two-sided
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$SSC,subset(nbq, Bact=="Xenic")$SSC, alternative = 'two.sided', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$SSC and subset(nbq, Bact == "Xenic")$SSC
## D = 0.09057, p-value = 0.0004998
## alternative hypothesis: two-sided
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$RED,subset(nbq, Bact=="Xenic")$RED, alternative = 'two.sided', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$RED and subset(nbq, Bact == "Xenic")$RED
## D = 0.49622, p-value = 0.0004998
## alternative hypothesis: two-sided
```

``` r
#NbO-D4
nbq=subset(read, strain=="NbO-D4")

ks.test(subset(nbq, Bact=="Axenic")$FSC,subset(nbq, Bact=="Xenic")$FSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$FSC and subset(nbq, Bact == "Xenic")$FSC
## D^+ = 0.069753, p-value = 0.0004998
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$SSC,subset(nbq, Bact=="Xenic")$SSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$SSC and subset(nbq, Bact == "Xenic")$SSC
## D^+ = 0.1177, p-value = 0.0004998
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$RED,subset(nbq, Bact=="Xenic")$RED, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$RED and subset(nbq, Bact == "Xenic")$RED
## D^+ = 0.23379, p-value = 0.0004998
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
#NbO-A5
nbq=subset(read, strain=="NbO-A5")

ks.test(subset(nbq, Bact=="Axenic")$FSC,subset(nbq, Bact=="Xenic")$FSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$FSC and subset(nbq, Bact == "Xenic")$FSC
## D^+ = 0.0020021, p-value = 0.9825
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$SSC,subset(nbq, Bact=="Xenic")$SSC, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$SSC and subset(nbq, Bact == "Xenic")$SSC
## D^+ = 0.0036998, p-value = 0.9485
## alternative hypothesis: the CDF of x lies above that of y
```

``` r
ks.test(subset(nbq, Bact=="Axenic")$RED,subset(nbq, Bact=="Xenic")$RED, alternative = 'greater', simulate.p.value=TRUE)
```

```
## 
## 	Monte-Carlo two-sample Kolmogorov-Smirnov test
## 
## data:  subset(nbq, Bact == "Axenic")$RED and subset(nbq, Bact == "Xenic")$RED
## D^+ = 0.0048131, p-value = 0.9115
## alternative hypothesis: the CDF of x lies above that of y
```
              
# 16S data

```
## Warning: package 'DESeq2' was built under R version 4.3.3
```

```
## Warning: package 'GenomeInfoDb' was built under R version 4.3.3
```

```
## Warning: package 'matrixStats' was built under R version 4.3.3
```

```
## Warning: package 'zCompositions' was built under R version 4.3.3
```

```
## Warning: package 'car' was built under R version 4.3.3
```

```
## [1] 1175   50
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1175 taxa and 50 samples ]
## sample_data() Sample Data:       [ 50 samples by 25 sample variables ]
## tax_table()   Taxonomy Table:    [ 1175 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1175 tips and 1174 internal nodes ]
```

```
##   OMA1  OMA11  OMA12  OMA13  OMA17  OMA18  OMA19   OMA2  OMA20  OMA21  OMA22 
##  30093  14214  16416   5655  24348   7968  19545  28086  26875  22995  29126 
##  OMA26  OMA27  OMA28   OMA3   OMA4  OMA47  OMA48  WGA10  WGA12  WGA13  WGA15 
##  12110   4862  17997  24888  31377  49845  23398  54055  15988  15713     19 
##  WGA18  WGA19  WGA20  WGA23  WGA25  WGA26  WGA27  WGA29   WGA3  WGA32  WGA33 
##  39669  16347   8355  15150  15915   1077  31290  24202 135475    783   8223 
##  WGA34  WGA35  WGA39   WGA4  WGA40  WGA41  WGA42  WGA45  WGA46  WGA48   WGA5 
##  19349   7798  29340  49137   5791  21777  35739  15876  34254  14517  41708 
##  WGA50  WGA53  WGA54  WGA55  WGA56   WGA7 
##  29196  16148  20534  17032     29    201
```

```
## [1] TRUE
```


# General Alpha statistics 


### figure 1 - no control 

``` r
cg_filt2=subset_samples(cg_filt, definition !="bacterioplankton")
cg_filt3=subset_samples(cg_filt2, diatom_control !="control")
cg_filt2=cg_filt3
s=specnumber(t(otu_table(cg_filt2)))


boxplot((s)~sample_data(cg_filt2)$treatment_26.1, las=2, xlab="", log='y')
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

``` r
treat2=sample_data(cg_filt2)$treatment_26.1

spec=cbind(s, treat2)
spec2=data.frame(spec)

spec2
```

```
##         s                   treat2
## OMA11 191             C1 Coculture
## OMA12 164             C1 Coculture
## OMA13 159             C1 Coculture
## OMA17 159             A1 Coculture
## OMA18 160             A1 Coculture
## OMA19 141             A1 Coculture
## OMA20 163             A2 Coculture
## OMA21 160             A2 Coculture
## OMA22 159             A2 Coculture
## OMA26 184             B1 Coculture
## OMA27 180             B1 Coculture
## OMA28 211             B1 Coculture
## WGA10  36 A1 Algal-cell associated
## WGA12  34 A1 Algal-cell associated
## WGA13  45 A1 Algal-cell associated
## WGA18  18 A2 Algal-cell associated
## WGA19  19 A2 Algal-cell associated
## WGA20  22 A2 Algal-cell associated
## WGA25  22 A2 Algal-cell associated
## WGA27  31 A2 Algal-cell associated
## WGA3   56 A1 Algal-cell associated
## WGA33  17 C1 Algal-cell associated
## WGA34  14 C1 Algal-cell associated
## WGA39  19 C1 Algal-cell associated
## WGA4   56 A1 Algal-cell associated
## WGA40  22 C1 Algal-cell associated
## WGA41  19 C1 Algal-cell associated
## WGA45  25 B1 Algal-cell associated
## WGA46  21 B1 Algal-cell associated
## WGA48  23 B1 Algal-cell associated
## WGA5   62 A1 Algal-cell associated
## WGA53  21 B1 Algal-cell associated
## WGA54  21 B1 Algal-cell associated
## WGA55  21 B1 Algal-cell associated
```

``` r
#mean1=aggregate(as.numeric(s) ~ treat3, data = spec2, FUN=mean)
mean1=with(spec2, tapply(as.numeric(s), treat2, mean))
sd1=with(spec2, tapply(as.numeric(s), treat2, sd))

treatment=names(sd1)
spec_bar=cbind(mean1, sd1)
spec_bar=data.frame(spec_bar)
spec_bar$treatment=treatment

custom_order =c( "A1 Coculture","A1 Algal-cell associated", 
                 "A2 Coculture","A2 Algal-cell associated",
                 "B1 Coculture","B1 Algal-cell associated",
                 "C1 Coculture","C1 Algal-cell associated")  
spec_bar$treatment2=spec_bar$treatment
spec_bar$treatment2 <- factor(spec_bar$treatment2, levels = custom_order)
spec_bar$treatment2=spec_bar$treatment
spec_bar$treatment2 <- factor(spec_bar$treatment2, levels = custom_order)

spec_bar$treatment_ag=spec_bar$treatment
spec_bar$treatment_ag=as.factor((spec_bar$treatment_ag))
levels(spec_bar$treatment_ag)=c("Algal-cell associated","Coculture", "Algal-cell associated","Coculture", "Algal-cell associated","Coculture", "Algal-cell associated", "Coculture")
spec_bar <- spec_bar %>%
  mutate(treatment_ag = fct_rev(factor(treatment_ag)))

spec_bar$name=spec_bar$treatment
spec_bar$name=as.factor((spec_bar$name))
levels(spec_bar$name)=c("A1", "A1", "A2", "A2", "B1", "B1", "C1", "C1")

values_colb=c( "#4D9DE0", '#E15554', '#3BB273', "#7768AE", "#4D9DE0", '#E15554', '#3BB273', "#7768AE")

values_colb=c( "#E64B35", "#4DBBD5", "#00A087", "#3C5488","#E64B35", "#4DBBD5", "#00A087", "#3C5488")


values_colb=c("#631879", "#3B4992", "#EE0000", "#008B45", "#631879", "#3B4992", "#EE0000", "#008B45")


## log which i like 
ggplot(spec_bar,aes(x=name, y=mean1)) +
  geom_bar( aes(x=name, y=mean1), stat="identity", fill=values_colb, alpha=0.7) +
    geom_errorbar( aes(x=name, ymin=mean1-sd1, ymax=mean1+sd1), width=0.4, colour="black", alpha=0.9, size=.7) +
  xlab(" ") +
  ylab("Species Richness") +
#  coord_flip() +
  facet_wrap(~treatment_ag, scales='free_x') +
  theme_classic(base_size=14) +
 scale_y_continuous(trans = 'log10') +
    theme(legend.position="bottom",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
```

```
## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
## ℹ Please use `linewidth` instead.
## This warning is displayed once per session.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-7-2.png)<!-- -->

``` r
figure1=ggplot(spec_bar,aes(x=name, y=mean1)) +
  geom_bar( aes(x=name, y=mean1), stat="identity", fill=values_colb, alpha=0.7) +
    geom_errorbar( aes(x=name, ymin=mean1-sd1, ymax=mean1+sd1), width=0.4, colour="black", alpha=0.9, size=.7) +
  xlab(" ") +
  ylab("Species Richness") +
#  coord_flip() +
  facet_wrap(~treatment_ag, scales='free_x') +
  theme_classic(base_size=14) +
 scale_y_continuous(trans = 'log10') +
    theme(legend.position="bottom",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=12)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#ggsave('figure1.svg',figure1, dpi=300, height=4, width=5)
```

# Figure 2 barplot

``` r
cg_filter=subset_samples(cg_filt,diatom_control=="diatom")
glom=tax_glom(cg_filter, "species")

abund=transform_sample_counts(glom, function(x) (x/sum(x)))

ps_filtered=ps_prune(abund, min.samples = 0, min.reads = 0, min.abundance = 0.015)
```

```
## 414 features grouped as 'Others' in the output
```

``` r
# 3. Melt the filtered data
df_melt <- psmelt(ps_filtered)
summary(as.factor(df_melt$species))
```

```
##         Aestuariibacter              Aquibacter Candidatus Pelagibacter 
##                      34                      34                      34 
##            Cellulophaga         Marinobacterium             Marinomonas 
##                      34                      34                      34 
##          Neptuniibacter             Owenweeksia           Tenacibaculum 
##                      34                      34                      34 
##                  Vibrio                Wandonia                    NA's 
##                      34                      34                      34
```

``` r
sub1=subset(df_melt, Sample == "WGA18")
sub1[,26:29]
```

```
##     keep_2026 diatom_control innoc  Kingdom
## 296       yes         diatom  <NA> Bacteria
## 84        yes         diatom  <NA> Bacteria
## 233       yes         diatom  <NA> Bacteria
## 19        yes         diatom  <NA>     <NA>
## 200       yes         diatom  <NA> Bacteria
## 360       yes         diatom  <NA> Bacteria
## 404       yes         diatom  <NA> Bacteria
## 329       yes         diatom  <NA> Bacteria
## 56        yes         diatom  <NA> Bacteria
## 128       yes         diatom  <NA> Bacteria
## 147       yes         diatom  <NA> Bacteria
## 269       yes         diatom  <NA> Bacteria
```

``` r
colors=c(
'Aestuariibacter' ="#e59d97",
'Aquibacter' ="#ffcaa5",
'Candidatus Pelagibacter' ="#e9daa2",
'Cellulophaga' ="#afb780",
'Marinobacterium' ="#f8fff1",
'Marinomonas' ="#8fc49b",
'Neptuniibacter' ="#93efe7",
'Owenweeksia' ="#ada5d9",
'Tenacibaculum' ="#f5e4ff",
'Vibrio' ="#ffc9e2",
'Wandonia' ="#c6a2ae",
'NA' ='white')


df_melt <- df_melt %>%
  mutate(definition_26 = fct_rev(factor(definition_26)))

#phy_df <- psmelt(abund)
p <- ggplot(df_melt, aes(x = Sample, y = Abundance, fill = species)) +
  theme_classic(base_size = 12) +
  xlab(" ") +
  ylab(" Relative Abundance") +
    facet_grid(~strain_26+definition_26,scales = "free", space="free") +
  geom_bar(stat = "identity",linewidth=0.3, color='black') +
  scale_fill_manual(values=colors) +
theme(axis.text.x = element_text(angle=90, size=10), legend.position = 'left') 
p
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


# nst

``` r
library(NST)

cg_filt2=subset_samples(cg_filt, definition !="bacterioplankton")

cg_filt2
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1028 taxa and 40 samples ]
## sample_data() Sample Data:       [ 40 samples by 25 sample variables ]
## tax_table()   Taxonomy Table:    [ 1028 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1028 tips and 1027 internal nodes ]
```

``` r
try1=subset_samples(cg_filt2, Treatment3!="Bact_T0_filtered")
try1=subset_samples(try1, Treatment3!="Bact_T0")
try1=subset_samples(try1, treatment_26!="Bacteria Control T5")

meta=cbind( sample_data(try1)$Treatment)
meta=data.frame(meta)
row.names(meta)=sample_data(try1)$X.OTU.ID
colnames(meta)=c( "Strain") # metagroup - strain and group - BA_Micro

groupie=cbind( sample_data(try1)$treatment_26)
groupie=data.frame(groupie)
row.names(groupie)=row.names(sample_data(try1))
colnames(groupie)=c( "Treatment") # metagroup - strain and group - BA_Micro
dim(groupie)
```

```
## [1] 34  1
```

``` r
otu_tab=t(otu_table(try1))
# taxonomic nst 
nst_1=tNST(comm = otu_tab , group = groupie, abundance.weighted=FALSE, meta.group = meta, dist.method = 'jaccard')
```

```
## All match very well.
```

```
## Warning in groupck(group): some groups have less than 6 samples, for which NST
## can be calculated but not recommened.
```

```
## Now randomizing by parallel computing. Begin at Wed Jul  1 12:28:37 2026. Please wait...
```

``` r
nst_data=nst_1$index.grp
par(mar=c(10,5,1,1))
groups_ag=c('cocultured','cocultured','cocultured','cocultured',"algal-cell associated","algal-cell associated","algal-cell associated","algal-cell associated")

groups_ag=c('co','co','co','co',"algal-cell","algal-cell","algal-cell","algal-cell")
boxplot(nst_data$NST.i.jaccard~groups_ag, xlab=" ", ylab="NST", ylim=c(0,1),yaxt='n')
axis(2, at=c(0,0.25,0.5,0.75,1))
box(which='plot')
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

``` r
nst_ba=nst_data$NST.i.jaccard[1:4]
nst_microbiome=nst_data$NST.i.jaccard[5:8]



mean(nst_ba)
```

```
## [1] 0.8480948
```

``` r
sd(nst_ba)
```

```
## [1] 0.03066064
```

``` r
mean(nst_microbiome)
```

```
## [1] 0.2099134
```

``` r
sd(nst_microbiome)
```

```
## [1] 0.06923988
```





``` r
library(NST)

cg_filt2=subset_samples(cg_filt, definition !="bacterioplankton")

cg_filt2
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1028 taxa and 40 samples ]
## sample_data() Sample Data:       [ 40 samples by 25 sample variables ]
## tax_table()   Taxonomy Table:    [ 1028 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1028 tips and 1027 internal nodes ]
```

``` r
try1=subset_samples(cg_filt2, Treatment3!="Bact_T0_filtered")
try1=subset_samples(try1, Treatment3!="Bact_T0")
try1=subset_samples(try1, treatment_26!="Bacteria Control T5")

strain=cbind( sample_data(try1)$diatom_strain)
strain=data.frame(strain)
row.names(strain)=sample_data(try1)$X.OTU.ID
colnames(strain)=c( "Strain") # metagroup - strain and group - BA_Micro

groupie=cbind( sample_data(try1)$definition)
groupie=data.frame(groupie)
row.names(groupie)=sample_data(try1)$X.OTU.ID
colnames(groupie)=c( "Treatment") # metagroup - strain and group - BA_Micro
dim(groupie)
```

```
## [1] 34  1
```

``` r
otu_tab=t(otu_table(try1))
# taxonomic nst 
nst_1=tNST(comm = otu_tab , group = groupie, meta.group = strain, dist.method = 'jaccard',abundance.weighted=FALSE)
```

```
## All match very well.
```

```
## Now randomizing by parallel computing. Begin at Wed Jul  1 12:28:47 2026. Please wait...
```

``` r
nst_data_ag=nst_1$index.grp
nst_data_ag
```

```
##        group size ST.i.jaccard NST.i.jaccard MST.i.jaccard
## 1 assemblage   66    0.8313426     0.7593887     0.7349763
## 2 microbiome  231    0.3063246     0.2888583     0.2610735
```

``` r
nst_data_group=nst_data_ag$ST.i.jaccard
names(nst_data_group)=c("co", "algal-cell")
par(mar=c(10,5,1,1))
groups_ag=c('co',"algal-cell")
barplot(nst_data_group,xlab=" ", ylab="NST", ylim=c(0,1),yaxt='n')
axis(2, at=c(0,0.25,0.5,0.75,1))
box(which='plot')
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


``` r
library(NST)

cg_filt2=subset_samples(cg_filt, definition !="bacterioplankton")

cg_filt2
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1028 taxa and 40 samples ]
## sample_data() Sample Data:       [ 40 samples by 25 sample variables ]
## tax_table()   Taxonomy Table:    [ 1028 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1028 tips and 1027 internal nodes ]
```

``` r
try1=subset_samples(cg_filt2, Treatment3!="Bact_T0_filtered")
try1=subset_samples(try1, Treatment3!="Bact_T0")
try1=subset_samples(try1, treatment_26!="Bacteria Control T5")

strain=cbind( sample_data(try1)$diatom_strain)
strain=data.frame(strain)
row.names(strain)=sample_data(try1)$X.OTU.ID
colnames(strain)=c( "Strain") # metagroup - strain and group - BA_Micro

groupie=cbind( sample_data(try1)$definition)
groupie=data.frame(groupie)
row.names(groupie)=sample_data(try1)$X.OTU.ID
colnames(groupie)=c( "Treatment") # metagroup - strain and group - BA_Micro
dim(groupie)
```

```
## [1] 34  1
```

``` r
otu_tab=t(otu_table(try1))
# taxonomic nst 
nst_1=tNST(comm = otu_tab , group = groupie, meta.group = strain, dist.method = 'jaccard',abundance.weighted=FALSE)
```

```
## All match very well.
```

```
## Now randomizing by parallel computing. Begin at Wed Jul  1 12:28:56 2026. Please wait...
```

``` r
nst_data_ag=nst_1$index.grp
nst_data_ag
```

```
##        group size ST.i.jaccard NST.i.jaccard MST.i.jaccard
## 1 assemblage   66    0.8309966     0.7616481     0.7346390
## 2 microbiome  231    0.3050131     0.2868303     0.2598057
```

``` r
nst_data_group=nst_data_ag$ST.i.jaccard
names(nst_data_group)=c("co", "algal-cell")
par(mar=c(7,5,1,1))
groups_ag=c('co',"algal-cell")
barplot(nst_data_group,xlab=" ", ylab="NST", ylim=c(0,1),yaxt='n', col=c("gray60","gray40"),las=2,cex.names=1.2)
axis(2, at=c(0,0.25,0.5,0.75,1))
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

``` r
#box(which='plot')
```
# Beta diversity and ordinations

## all samples pcoa

``` r
diatom=subset_samples(cg_filt, diatom_control=='diatom')

diatom2=transform_sample_counts(diatom, function(x) log(x+1))
bray=vegdist(t(otu_table(diatom2)),'bray')

p=prcomp(bray)
summary(p)
```

```
## Importance of components:
##                           PC1    PC2     PC3     PC4     PC5     PC6     PC7
## Standard deviation     0.7146 0.3310 0.21600 0.19124 0.15878 0.15025 0.13917
## Proportion of Variance 0.5480 0.1176 0.05008 0.03925 0.02706 0.02423 0.02079
## Cumulative Proportion  0.5480 0.6656 0.71571 0.75497 0.78203 0.80626 0.82704
##                            PC8     PC9    PC10    PC11    PC12    PC13    PC14
## Standard deviation     0.13379 0.12590 0.11265 0.10715 0.10498 0.09774 0.09653
## Proportion of Variance 0.01921 0.01701 0.01362 0.01232 0.01183 0.01025 0.01000
## Cumulative Proportion  0.84626 0.86327 0.87689 0.88921 0.90104 0.91130 0.92130
##                           PC15    PC16    PC17    PC18    PC19    PC20    PC21
## Standard deviation     0.08530 0.08205 0.07601 0.07209 0.07058 0.06944 0.06744
## Proportion of Variance 0.00781 0.00723 0.00620 0.00558 0.00535 0.00518 0.00488
## Cumulative Proportion  0.92911 0.93633 0.94253 0.94811 0.95346 0.95864 0.96352
##                           PC22    PC23    PC24    PC25    PC26    PC27    PC28
## Standard deviation     0.06343 0.06287 0.06084 0.05944 0.05481 0.05256 0.05135
## Proportion of Variance 0.00432 0.00424 0.00397 0.00379 0.00322 0.00297 0.00283
## Cumulative Proportion  0.96783 0.97208 0.97605 0.97984 0.98307 0.98603 0.98886
##                           PC29    PC30    PC31    PC32    PC33      PC34
## Standard deviation     0.04973 0.04878 0.04574 0.04316 0.03961 4.468e-17
## Proportion of Variance 0.00265 0.00255 0.00225 0.00200 0.00168 0.000e+00
## Cumulative Proportion  0.99152 0.99407 0.99632 0.99832 1.00000 1.000e+00
```

``` r
{
  
  par(mar=c(10,
            5,1,1),xpd=FALSE)
  plot(p$x[,1],p$x[,2],pch=sample_data(diatom)$pch,
     bg=as.character(sample_data(diatom)$col_26),
     xlab = "PCoA1 54 %", 
     ylab= "PCoA2 12 %")
 # axis(2, at=c(-2,0,2))
ordiellipse(p$x,          group=as.factor(sample_data(diatom)$treatment_26.1),

           kind ='sd',conf=0.8,
           label=TRUE)
legend(-2, -5,
       legend=c("NbO-A5 BA","NbO-A5 Micro", "NbO-D4 BA", "NbO-D4 Micro", "NbP-YE5 BA", "NbP-YE5 Micro", "NbQ-B7 BA",
                "NbQ-B7 Micro"),
       bty='n',
       pt.bg=c("#4D9DE0", "#4D9DE0", "#E15554", "#E15554", "#3BB273", "#3BB273", "#7768AE", "#7768AE"),

       pch=c(21,22,21,22,21,22,21,22))

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/beta-1-1.png)<!-- -->

## just microbiomes

``` r
cg_filt2=subset_samples(cg_filt, definition_26 =='Algal-cell associated')
diatom_micro=transform_sample_counts(cg_filt2, function(x) log(x+1))
bray=vegdist(t(otu_table(diatom_micro)),'bray')

p_micro=prcomp(bray)
summary(p_micro)
```

```
## Importance of components:
##                           PC1    PC2    PC3    PC4     PC5     PC6     PC7
## Standard deviation     0.4481 0.2714 0.2498 0.1987 0.19204 0.17183 0.15999
## Proportion of Variance 0.3256 0.1194 0.1012 0.0640 0.05978 0.04787 0.04149
## Cumulative Proportion  0.3256 0.4450 0.5461 0.6101 0.66992 0.71778 0.75928
##                            PC8     PC9    PC10    PC11    PC12    PC13    PC14
## Standard deviation     0.15547 0.14003 0.12853 0.12094 0.10744 0.10379 0.09115
## Proportion of Variance 0.03919 0.03179 0.02678 0.02371 0.01871 0.01746 0.01347
## Cumulative Proportion  0.79846 0.83025 0.85703 0.88074 0.89945 0.91691 0.93038
##                           PC15    PC16    PC17    PC18    PC19    PC20    PC21
## Standard deviation     0.08961 0.08421 0.08363 0.07950 0.07519 0.06819 0.06486
## Proportion of Variance 0.01302 0.01150 0.01134 0.01024 0.00916 0.00754 0.00682
## Cumulative Proportion  0.94340 0.95490 0.96623 0.97648 0.98564 0.99318 1.00000
##                             PC22
## Standard deviation     6.306e-17
## Proportion of Variance 0.000e+00
## Cumulative Proportion  1.000e+00
```

``` r
mod <- betadisper(bray,group=as.factor(sample_data(diatom_micro)$strain_26))
anova(mod)
```

```
## Analysis of Variance Table
## 
## Response: Distances
##           Df   Sum Sq   Mean Sq F value Pr(>F)  
## Groups     3 0.042385 0.0141282  4.2966 0.0188 *
## Residuals 18 0.059188 0.0032882                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
permutest(mod)
```

```
## 
## Permutation test for homogeneity of multivariate dispersions
## Permutation: free
## Number of permutations: 999
## 
## Response: Distances
##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
## Groups     3 0.042385 0.0141282 4.2966    999  0.019 *
## Residuals 18 0.059188 0.0032882                       
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
TukeyHSD(mod)
```

```
##   Tukey multiple comparisons of means
##     95% family-wise confidence level
## 
## Fit: aov(formula = distances ~ group, data = df)
## 
## $group
##               diff         lwr       upr     p adj
## A2-A1  0.064367286 -0.03377008 0.1625046 0.2819977
## B1-A1  0.061581725 -0.03198858 0.1551520 0.2793144
## C1-A1  0.124306213  0.02616885 0.2224436 0.0104861
## B1-A2 -0.002785561 -0.10092292 0.0953518 0.9998075
## C1-A2  0.059938926 -0.04256221 0.1624401 0.3760067
## C1-B1  0.062724488 -0.03541287 0.1608618 0.3025753
```

``` r
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/beta-2a-1.png)<!-- -->

``` r
color_micro=c(sample_data(diatom_micro)$col_26)

{
  par(mar=c(10,
            5,1,1),xpd=FALSE)
  plot(p_micro$x[,1],p_micro$x[,2],pch=sample_data(diatom_micro)$pch,
     bg=as.character(sample_data(diatom_micro)$col_26), cex=2,
     xlab = "PCoA1 32.56 %", 
     ylab= "PCoA2 11.94 %")
 # axis(2, at=c(-2,0,2))
ordiellipse(p_micro$x,          group=as.factor(sample_data(diatom_micro)$strain_26),

           kind ='sd',conf=0.8,
           label=TRUE)
legend(-1, 0,
       legend=c("A1", "A2", "B1", "C1"),
       bty='n',
       pt.bg=c("#7768AE","#4D9DE0", "#E15554",  "#3BB273" ),

       pch=c(21,22,21,22,21,22,21,22))

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/beta-2a-2.png)<!-- -->

``` r
anosim(bray, sample_data(diatom_micro)$strain_26)
```

```
## 
## Call:
## anosim(x = bray, grouping = sample_data(diatom_micro)$strain_26) 
## Dissimilarity: bray 
## 
## ANOSIM statistic R: 0.6365 
##       Significance: 0.001 
## 
## Permutation: free
## Number of permutations: 999
```

``` r
adonis2(bray~sample_data(diatom_micro)$strain_26)
```

```
## Permutation test for adonis under reduced model
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = bray ~ sample_data(diatom_micro)$strain_26)
##          Df SumOfSqs      R2      F Pr(>F)    
## Model     3   1.6684 0.36128 3.3938  0.001 ***
## Residual 18   2.9496 0.63872                  
## Total    21   4.6180 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


## just cos

``` r
cg_filt2=subset_samples(cg_filt, definition_26 =='Coculture')
diatom_co=transform_sample_counts(cg_filt2, function(x) log(x+1))
bray=vegdist(t(otu_table(diatom_co)),'bray')

p_co=prcomp(bray)
summary(p_co)
```

```
## Importance of components:
##                           PC1    PC2    PC3     PC4     PC5     PC6     PC7
## Standard deviation     0.2276 0.1711 0.1255 0.11889 0.10708 0.10007 0.09486
## Proportion of Variance 0.3091 0.1747 0.0940 0.08438 0.06844 0.05978 0.05371
## Cumulative Proportion  0.3091 0.4839 0.5779 0.66224 0.73068 0.79046 0.84418
##                            PC8     PC9    PC10    PC11      PC12
## Standard deviation     0.09047 0.08259 0.07699 0.07191 6.193e-17
## Proportion of Variance 0.04886 0.04072 0.03538 0.03087 0.000e+00
## Cumulative Proportion  0.89304 0.93375 0.96913 1.00000 1.000e+00
```

``` r
mod <- betadisper(bray,group=as.factor(sample_data(diatom_co)$strain_26))
anova(mod)
```

```
## Analysis of Variance Table
## 
## Response: Distances
##           Df    Sum Sq    Mean Sq F value Pr(>F)
## Groups     3 0.0022242 0.00074139  0.3137 0.8153
## Residuals  8 0.0189093 0.00236367
```

``` r
permutest(mod)
```

```
## 
## Permutation test for homogeneity of multivariate dispersions
## Permutation: free
## Number of permutations: 999
## 
## Response: Distances
##           Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
## Groups     3 0.0022242 0.00074139 0.3137    999  0.798
## Residuals  8 0.0189093 0.00236367
```

``` r
TukeyHSD(mod)
```

```
##   Tukey multiple comparisons of means
##     95% family-wise confidence level
## 
## Fit: aov(formula = distances ~ group, data = df)
## 
## $group
##                diff         lwr       upr     p adj
## A2-A1  0.0326619034 -0.09445885 0.1597827 0.8423184
## B1-A1  0.0330131620 -0.09410759 0.1601339 0.8382356
## C1-A1  0.0273634830 -0.09975727 0.1544842 0.8983450
## B1-A2  0.0003512586 -0.12676949 0.1274720 0.9999997
## C1-A2 -0.0052984205 -0.13241917 0.1218223 0.9990771
## C1-B1 -0.0056496791 -0.13277043 0.1214711 0.9988826
```

``` r
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/beta-2-1.png)<!-- -->

``` r
{
  par(mar=c(10,
            5,1,1),xpd=FALSE)
  plot(p_co$x[,1],p_co$x[,2],pch=sample_data(diatom_co)$pch,
     bg=sample_data(diatom_co)$col_26, cex=2,
     xlab = "PCoA1 30.91 %", 
     ylab= "PCoA2 17.47 %")
 # axis(2, at=c(-2,0,2))
ordiellipse(p_co$x,          group=as.factor(sample_data(diatom_co)$strain_26),

           kind ='sd',conf=0.8,
           label=TRUE)
legend(-1, 0,
       legend=c("A1", "A2", "B1", "C1"),
       bty='n',
       pt.bg=c("#7768AE","#4D9DE0", "#E15554",  "#3BB273" ),

       pch=c(21,22,21,22,21,22,21,22))

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/beta-2-2.png)<!-- -->


## both pcoas together 

``` r
{
  par(mar=c(10,
            5,1,1),mfrow=c(1,2),xpd=FALSE)
  
  
  plot(p_co$x[,1],p_co$x[,2],pch=sample_data(diatom_co)$pch,
     bg=sample_data(diatom_co)$col_26, cex=2,
     xlab = "PCoA1 30.91 %", 
     ylab= "PCoA2 17.47 %", xaxt='none', yaxt='none')
  axis(1, at=c(-0.3,0,0.3))
  axis(2, at=c(-0.2,0,0.2))
 # axis(2, at=c(-2,0,2))
ordiellipse(p_co$x,          group=as.factor(sample_data(diatom_co)$strain_26),

           kind ='sd',conf=0.8,
           label=TRUE)
legend(-1, 0,
       legend=c("A1", "A2", "B1", "C1"),
       bty='n',
       pt.bg=c("#7768AE","#4D9DE0", "#E15554",  "#3BB273" ),

       pch=c(21,22,21,22,21,22,21,22))

plot(p_micro$x[,1],p_micro$x[,2],pch=sample_data(diatom_micro)$pch,
     bg=as.character(sample_data(diatom_micro)$col_26), cex=2,
     xlab = "PCoA1 32.56 %", 
     ylab= "PCoA2 11.94 %", xaxt='none', yaxt='none')
  axis(1, at=c(-0.5,0,0.5))
  axis(2, at=c(-0.5,0,0.5))
 # axis(2, at=c(-2,0,2))
ordiellipse(p_micro$x,          group=as.factor(sample_data(diatom_micro)$strain_26),

           kind ='sd',conf=0.8,
           label=TRUE)

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/pcoa-3-1.png)<!-- -->

``` r
  par(mar=c(20,5,1,1),mfrow=c(1,2),xpd=TRUE)
  plot(p_co$x[,1],p_co$x[,2],pch=sample_data(diatom_co)$pch,
     bg=sample_data(diatom_co)$col_26, cex=2,
     xlab = "PCoA1 30.91 %", 
     ylab= "PCoA2 17.47 %", xaxt='none', yaxt='none')
  axis(1, at=c(-0.3,0,0.3))
  axis(2, at=c(-0.2,0,0.2))
 # axis(2, at=c(-2,0,2))
legend(0, -0.7,
       legend=c("A1", "A2", "B1", "C1"),
       bty='n',
       pt.bg=c("#7768AE","#4D9DE0", "#E15554",  "#3BB273" ), pch=21,
       title='strain')
legend(0, -1.5,
       legend=c("Coculture", "Algal-cell Associated"),
       bty='n',
       pt.bg=c('white' ), pch=c(21,23),
       title='treatment')
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/pcoa-3-2.png)<!-- -->


### ggplot 2
#### microbiomes 


``` r
library(vegan)
library(ggplot2)
library(dplyr)

# ---- 1. Compute distance matrix ----
# otu_table: samples as rows, taxa/ASVs as columns (relative abundance or counts)
cg_filt2=subset_samples(cg_filt,  definition_26 =='Algal-cell associated')
diatom_co=transform_sample_counts(cg_filt2, function(x) log(x+1))
bray=vegdist(t(otu_table(diatom_co)),'bray')


# ---- 2. Run PCoA ----
pcoa_result <- prcomp(bray)
# Extract axis scores
pca_scores <- as.data.frame(pcoa_result$x[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$SampleID <- rownames(pca_scores)

# % variance explained per axis
var_explained <- round((pcoa_result$sdev^2 / sum(pcoa_result$sdev^2)) * 100, 1)


strain_colors=c("A1"="#7768AE",
                "A2"="#4D9DE0", 
                "B1"="#E15554", "C1"= "#3BB273" )

# ---- 3. Merge with sample metadata ----
# metadata should have a SampleID column matching otu_table rownames,
# plus the variable(s) you want to facet/color by
metadata=sample_data(cg_filt2)
metadata$SampleID=metadata$X.OTU.ID
pca_scores <- pca_scores %>%
  left_join(metadata, by = "SampleID")



calc_ellipse <- function(df, npoints = 100) {
  if (nrow(df) < 3) return(NULL)  # true mathematical floor: need ≥3 points for a non-degenerate 2D covariance matrix
  center <- colMeans(df[, c("PC1", "PC2")])
  cov_mat <- cov(df[, c("PC1", "PC2")])
  
  theta <- (0:npoints) * 2 * pi / npoints
  circle <- cbind(cos(theta), sin(theta))
  
  chol_decomp <- tryCatch(chol(cov_mat), error = function(e) NULL)
  if (is.null(chol_decomp)) return(NULL)  # singular matrix (e.g. n=2, points collinear)
  
  ellipse_pts <- t(center + t(circle %*% chol_decomp) * sqrt(qchisq(0.95, df = 2)))
  data.frame(PC1 = ellipse_pts[, 1], PC2 = ellipse_pts[, 2])
}

ellipse_df <- pca_scores %>%
  group_by( strain_26) %>%
  group_modify(~ calc_ellipse(.x)) %>%
  ungroup()



# ---- 4. Plot ----
pcoa1=ggplot(pca_scores, aes(x = PC1, y = PC2, fill = strain_26)) +
  geom_point(size = 3, alpha = 0.8, pch=23) +
  labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "Group"
  ) +
  theme_classic2(base_size = 14) +
  geom_path(data = ellipse_df, aes(x = PC1, y = PC2),
            linewidth = 0.3) +
  scale_fill_manual(values=strain_colors) + theme(legend.position = 'none')
pcoa1
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

#### coculutres 




``` r
library(vegan)
library(ggplot2)
library(dplyr)

# ---- 1. Compute distance matrix ----
# otu_table: samples as rows, taxa/ASVs as columns (relative abundance or counts)
cg_filt2=subset_samples(cg_filt,  definition_26 =='Coculture')
diatom_co=transform_sample_counts(cg_filt2, function(x) log(x+1))
bray=vegdist(t(otu_table(diatom_co)),'bray')


# ---- 2. Run PCoA ----
pcoa_result <- prcomp(bray)
# Extract axis scores
pca_scores <- as.data.frame(pcoa_result$x[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$SampleID <- rownames(pca_scores)

# % variance explained per axis
var_explained <- round((pcoa_result$sdev^2 / sum(pcoa_result$sdev^2)) * 100, 1)


strain_colors=c("A1"="#7768AE",
                "A2"="#4D9DE0", 
                "B1"="#E15554", "C1"= "#3BB273" )

# ---- 3. Merge with sample metadata ----
# metadata should have a SampleID column matching otu_table rownames,
# plus the variable(s) you want to facet/color by
metadata=sample_data(cg_filt2)
metadata$SampleID=metadata$X.OTU.ID
pca_scores <- pca_scores %>%
  left_join(metadata, by = "SampleID")

calc_ellipse <- function(df, npoints = 100) {
  if (nrow(df) < 3) return(NULL)  # true mathematical floor: need ≥3 points for a non-degenerate 2D covariance matrix
  center <- colMeans(df[, c("PC1", "PC2")])
  cov_mat <- cov(df[, c("PC1", "PC2")])
  
  theta <- (0:npoints) * 2 * pi / npoints
  circle <- cbind(cos(theta), sin(theta))
  
  chol_decomp <- tryCatch(chol(cov_mat), error = function(e) NULL)
  if (is.null(chol_decomp)) return(NULL)  # singular matrix (e.g. n=2, points collinear)
  
  ellipse_pts <- t(center + t(circle %*% chol_decomp) * sqrt(qchisq(0.95, df = 2)))
  data.frame(PC1 = ellipse_pts[, 1], PC2 = ellipse_pts[, 2])
}

ellipse_df <- pca_scores %>%
  group_by( strain_26) %>%
  group_modify(~ calc_ellipse(.x)) %>%
  ungroup()

pcoa2=ggplot(pca_scores, aes(x = PC1, y = PC2, fill = strain_26)) +
  geom_point(size = 3, alpha = 0.8, pch=23) +
  labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "Group"
  ) +
  theme_classic2(base_size = 14) +
  geom_path(data = ellipse_df, aes(x = PC1, y = PC2),
            linewidth = 0.3) +
  scale_fill_manual(values=strain_colors) + theme(legend.position = 'none')

pcoa2
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-13-1.png)<!-- -->



``` r
cows=cowplot::plot_grid(pcoa2,pcoa1,labels='auto')
#ggsave('pcoa-ggplot.svg', cows, height=4, width=10)
```


#### legend 


``` r
library(vegan)
library(ggplot2)
library(dplyr)

# ---- 1. Compute distance matrix ----
# otu_table: samples as rows, taxa/ASVs as columns (relative abundance or counts)
cg_filt2=subset_samples(cg_filt,  definition_26 =='Coculture' |
                           definition_26 =='Algal-cell associated')
diatom_co=transform_sample_counts(cg_filt2, function(x) log(x+1))
bray=vegdist(t(otu_table(diatom_co)),'bray')


# ---- 2. Run PCoA ----
pcoa_result <- prcomp(bray)
# Extract axis scores
pca_scores <- as.data.frame(pcoa_result$x[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$SampleID <- rownames(pca_scores)

# % variance explained per axis
var_explained <- round((pcoa_result$sdev^2 / sum(pcoa_result$sdev^2)) * 100, 1)


strain_colors=c("A1"="#7768AE",
                "A2"="#4D9DE0", 
                "B1"="#E15554", "C1"= "#3BB273" )

# ---- 3. Merge with sample metadata ----
# metadata should have a SampleID column matching otu_table rownames,
# plus the variable(s) you want to facet/color by
metadata=sample_data(cg_filt2)
metadata$SampleID=metadata$X.OTU.ID
pca_scores <- pca_scores %>%
  left_join(metadata, by = "SampleID")

calc_ellipse <- function(df, npoints = 100) {
  if (nrow(df) < 3) return(NULL)  # true mathematical floor: need ≥3 points for a non-degenerate 2D covariance matrix
  center <- colMeans(df[, c("PC1", "PC2")])
  cov_mat <- cov(df[, c("PC1", "PC2")])
  
  theta <- (0:npoints) * 2 * pi / npoints
  circle <- cbind(cos(theta), sin(theta))
  
  chol_decomp <- tryCatch(chol(cov_mat), error = function(e) NULL)
  if (is.null(chol_decomp)) return(NULL)  # singular matrix (e.g. n=2, points collinear)
  
  ellipse_pts <- t(center + t(circle %*% chol_decomp) * sqrt(qchisq(0.95, df = 2)))
  data.frame(PC1 = ellipse_pts[, 1], PC2 = ellipse_pts[, 2])
}

ellipse_df <- pca_scores %>%
  group_by( strain_26) %>%
  group_modify(~ calc_ellipse(.x)) %>%
  ungroup()


pca_scores_all=pca_scores


pcoa_legs=ggplot(pca_scores_all, aes(x = PC1, y = PC2, fill = strain_26, color=strain_26)) +
  geom_point(size = 3, alpha = 0.8,aes(fill=strain_26, shape=definition_26)) +
  labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "Group"
  ) +
  facet_wrap(~definition_26) +
  theme_classic2(base_size = 14) +
  scale_fill_manual(values=strain_colors) +
  scale_color_manual(values=strain_colors)+ theme(legend.position = 'left')
#ggsave('pcoa-legend.svg', plot=pcoa_legs)
```


# ribbon plot


``` r
library(ggalluvial)
library(ggplot2)
df_melt
```

```
##        OTU Sample    Abundance X.OTU.ID      info strain_26 diatom_strain
## 296   Otu4  WGA18 0.6501439737    WGA18 WGA_RD4-A        A2        NbO-D4
## 295   Otu4  WGA19 0.6147675854    WGA19 WGA_RD4-A        A2        NbO-D4
## 281   Otu4  OMA22 0.5571781160    OMA22  RD4-C_18        A2        NbO-D4
## 111  Otu12  WGA33 0.5277663513    WGA33 WGA_NB7-A        C1        NbQ-B7
## 18  Others  WGA40 0.5168654875    WGA40 WGA_NB7-B        C1        NbQ-B7
## 224   Otu2  WGA34 0.5143973673    WGA34 WGA_NB7-A        C1        NbQ-B7
## 231   Otu2  OMA20 0.5139860140    OMA20  RD4-A_18        A2        NbO-D4
## 222   Otu2  OMA13 0.5111065405    OMA13  NB7-C_18        C1        NbQ-B7
## 7   Others  WGA39 0.5074043603    WGA39 WGA_NB7-B        C1        NbQ-B7
## 279   Otu4  WGA20 0.4845742493    WGA20 WGA_RD4-A        A2        NbO-D4
## 214   Otu2  OMA12 0.4790209790    OMA12  NB7-B_18        C1        NbQ-B7
## 358   Otu6  WGA41 0.4625668449    WGA41 WGA_NB7-B        C1        NbQ-B7
## 390   Otu7  OMA19 0.4574249280    OMA19  RA5-C_18        A1        NbO-A5
## 205   Otu2  OMA26 0.4477581242    OMA26  YE5-A_18        B1       NbP-YE5
## 207   Otu2  OMA27 0.4224598930    OMA27  YE5-B_18        B1       NbP-YE5
## 289   Otu4  WGA53 0.3879062114    WGA53 WGA_YE5-B        B1       NbP-YE5
## 362   Otu6  WGA54 0.3835870012    WGA54 WGA_YE5-B        B1       NbP-YE5
## 206   Otu2  OMA11 0.3745372275    OMA11  NB7-A_18        C1        NbQ-B7
## 406   Otu7   WGA4 0.3712464007     WGA4 WGA_RA5-A        A1        NbO-A5
## 299   Otu4  WGA10 0.3679555738    WGA10 WGA_RA5-B        A1        NbO-A5
## 12  Others  OMA28 0.3455368161    OMA28  YE5-C_18        B1       NbP-YE5
## 290   Otu4  WGA27 0.3443027561    WGA27 WGA_RD4-B        A2        NbO-D4
## 379   Otu7  OMA17 0.3414232826    OMA17  RA5-A_18        A1        NbO-A5
## 211   Otu2  OMA28 0.3364870424    OMA28  YE5-C_18        B1       NbP-YE5
## 34  Others  WGA55 0.3362813657    WGA55 WGA_YE5-B        B1       NbP-YE5
## 141  Otu14  WGA25 0.3354586590    WGA25 WGA_RD4-B        A2        NbO-D4
## 220   Otu2  OMA18 0.3280542986    OMA18  RA5-B_18        A1        NbO-A5
## 210   Otu2  OMA21 0.3276429453    OMA21  RD4-B_18        A2        NbO-D4
## 99   Otu11  WGA34 0.3175647882    WGA34 WGA_NB7-A        C1        NbQ-B7
## 326   Otu5  WGA46 0.3171534348    WGA46 WGA_YE5-A        B1       NbP-YE5
## 384   Otu7   WGA3 0.3140682847     WGA3 WGA_RA5-A        A1        NbO-A5
## 292   Otu4  WGA46 0.3138626080    WGA46 WGA_YE5-A        B1       NbP-YE5
## 5   Others  OMA11 0.3027560675    OMA11  NB7-A_18        C1        NbQ-B7
## 13  Others  OMA12 0.2982311806    OMA12  NB7-B_18        C1        NbQ-B7
## 256  Otu24  WGA48 0.2914438503    WGA48 WGA_YE5-A        B1       NbP-YE5
## 305   Otu4  WGA55 0.2815713698    WGA55 WGA_YE5-B        B1       NbP-YE5
## 181  Otu19  WGA45 0.2782805430    WGA45 WGA_YE5-A        B1       NbP-YE5
## 277   Otu4  OMA21 0.2747840395    OMA21  RD4-B_18        A2        NbO-D4
## 28  Others  OMA26 0.2556561086    OMA26  YE5-A_18        B1       NbP-YE5
## 1   Others   WGA5 0.2542163719     WGA5 WGA_RA5-A        A1        NbO-A5
## 2   Others  WGA45 0.2519539284    WGA45 WGA_YE5-A        B1       NbP-YE5
## 30  Others  WGA41 0.2441382147    WGA41 WGA_NB7-B        C1        NbQ-B7
## 32  Others  OMA20 0.2414644179    OMA20  RD4-A_18        A2        NbO-D4
## 168  Otu14  WGA41 0.2400246812    WGA41 WGA_NB7-B        C1        NbQ-B7
## 21  Others  OMA18 0.2398190045    OMA18  RA5-B_18        A1        NbO-A5
## 235   Otu2  WGA12 0.2381735911    WGA12 WGA_RA5-B        A1        NbO-A5
## 4   Others  OMA27 0.2354997943    OMA27  YE5-B_18        B1       NbP-YE5
## 10  Others  WGA53 0.2334430276    WGA53 WGA_YE5-B        B1       NbP-YE5
## 86   Otu11  WGA12 0.2305635541    WGA12 WGA_RA5-B        A1        NbO-A5
## 25  Others  OMA13 0.2295351707    OMA13  NB7-C_18        C1        NbQ-B7
## 20  Others  WGA13 0.2295351707    WGA13 WGA_RA5-B        A1        NbO-A5
## 387   Otu7  WGA25 0.2266556972    WGA25 WGA_RD4-B        A2        NbO-D4
## 208   Otu2  OMA17 0.2225421637    OMA17  RA5-A_18        A1        NbO-A5
## 218   Otu2  WGA27 0.2219251337    WGA27 WGA_RD4-B        A2        NbO-D4
## 236   Otu2  WGA39 0.2219251337    WGA39 WGA_NB7-B        C1        NbQ-B7
## 14  Others  WGA46 0.2186343069    WGA46 WGA_YE5-A        B1       NbP-YE5
## 26  Others  WGA48 0.2159605101    WGA48 WGA_YE5-A        B1       NbP-YE5
## 219   Otu2  OMA19 0.2108185932    OMA19  RA5-C_18        A1        NbO-A5
## 398   Otu7   WGA5 0.2104072398     WGA5 WGA_RA5-A        A1        NbO-A5
## 9   Others  OMA17 0.2001234060    OMA17  RA5-A_18        A1        NbO-A5
## 288   Otu4  WGA33 0.1995063760    WGA33 WGA_NB7-A        C1        NbQ-B7
## 24  Others  WGA10 0.1879884821    WGA10 WGA_RA5-B        A1        NbO-A5
## 276   Otu4   WGA3 0.1863430687     WGA3 WGA_RA5-A        A1        NbO-A5
## 228   Otu2  WGA48 0.1803784451    WGA48 WGA_YE5-A        B1       NbP-YE5
## 223   Otu2  WGA10 0.1793500617    WGA10 WGA_RA5-B        A1        NbO-A5
## 22  Others  WGA54 0.1762649116    WGA54 WGA_YE5-B        B1       NbP-YE5
## 3   Others   WGA3 0.1737967914     WGA3 WGA_RA5-A        A1        NbO-A5
## 199  Otu19  WGA19 0.1723570547    WGA19 WGA_RD4-A        A2        NbO-D4
## 332   Otu5  WGA12 0.1713286713    WGA12 WGA_RA5-B        A1        NbO-A5
## 17  Others  OMA21 0.1709173180    OMA21  RD4-B_18        A2        NbO-D4
## 16  Others  WGA20 0.1707116413    WGA20 WGA_RD4-A        A2        NbO-D4
## 216   Otu2  WGA20 0.1655697244    WGA20 WGA_RD4-A        A2        NbO-D4
## 337   Otu5   WGA5 0.1631016043     WGA5 WGA_RA5-A        A1        NbO-A5
## 230   Otu2  WGA40 0.1606334842    WGA40 WGA_NB7-B        C1        NbQ-B7
## 33  Others  OMA19 0.1585767174    OMA19  RA5-C_18        A1        NbO-A5
## 226   Otu2  WGA53 0.1559029206    WGA53 WGA_YE5-B        B1       NbP-YE5
## 203  Otu19  WGA13 0.1499382970    WGA13 WGA_RA5-B        A1        NbO-A5
## 354   Otu6  WGA48 0.1497326203    WGA48 WGA_YE5-A        B1       NbP-YE5
## 84   Otu11  WGA18 0.1495269436    WGA18 WGA_RD4-A        A2        NbO-D4
## 117  Otu12  WGA45 0.1472645002    WGA45 WGA_YE5-A        B1       NbP-YE5
## 215   Otu2  WGA25 0.1470588235    WGA25 WGA_RD4-B        A2        NbO-D4
## 27  Others  WGA27 0.1460304401    WGA27 WGA_RD4-B        A2        NbO-D4
## 293   Otu4  WGA13 0.1431509667    WGA13 WGA_RA5-B        A1        NbO-A5
## 98   Otu11  WGA40 0.1404771699    WGA40 WGA_NB7-B        C1        NbQ-B7
## 233   Otu2  WGA18 0.1402714932    WGA18 WGA_RD4-A        A2        NbO-D4
## 8   Others  WGA12 0.1394487865    WGA12 WGA_RA5-B        A1        NbO-A5
## 131  Otu12  WGA39 0.1375976964    WGA39 WGA_NB7-B        C1        NbQ-B7
## 306   Otu4   WGA4 0.1367749897     WGA4 WGA_RA5-A        A1        NbO-A5
## 186  Otu19  OMA18 0.1330728095    OMA18  RA5-B_18        A1        NbO-A5
## 29  Others  OMA22 0.1330728095    OMA22  RD4-C_18        A2        NbO-D4
## 213   Otu2   WGA3 0.1314273961     WGA3 WGA_RA5-A        A1        NbO-A5
## 134  Otu12  WGA53 0.1312217195    WGA53 WGA_YE5-B        B1       NbP-YE5
## 95   Otu11   WGA4 0.1287535993     WGA4 WGA_RA5-A        A1        NbO-A5
## 238   Otu2   WGA4 0.1256684492     WGA4 WGA_RA5-A        A1        NbO-A5
## 322   Otu5  WGA54 0.1248457425    WGA54 WGA_YE5-B        B1       NbP-YE5
## 393   Otu7  WGA12 0.1246400658    WGA12 WGA_RA5-B        A1        NbO-A5
## 391   Otu7  OMA18 0.1229946524    OMA18  RA5-B_18        A1        NbO-A5
## 343   Otu6  WGA45 0.1225832991    WGA45 WGA_YE5-A        B1       NbP-YE5
## 227   Otu2   WGA5 0.1223776224     WGA5 WGA_RA5-A        A1        NbO-A5
## 297   Otu4  WGA54 0.1209378856    WGA54 WGA_YE5-B        B1       NbP-YE5
## 15  Others  WGA25 0.1194981489    WGA25 WGA_RD4-B        A2        NbO-D4
## 148  Otu14  WGA13 0.1180584122    WGA13 WGA_RA5-B        A1        NbO-A5
## 237   Otu2  WGA55 0.1164129988    WGA55 WGA_YE5-B        B1       NbP-YE5
## 370   Otu6  WGA55 0.1129164953    WGA55 WGA_YE5-B        B1       NbP-YE5
## 217   Otu2  WGA45 0.1125051419    WGA45 WGA_YE5-A        B1       NbP-YE5
## 394   Otu7  WGA10 0.1120937886    WGA10 WGA_RA5-B        A1        NbO-A5
## 209   Otu2  OMA22 0.1116824352    OMA22  RD4-C_18        A2        NbO-D4
## 313   Otu5  WGA33 0.1104483752    WGA33 WGA_NB7-A        C1        NbQ-B7
## 232   Otu2  WGA19 0.1071575483    WGA19 WGA_RD4-A        A2        NbO-D4
## 174  Otu19  OMA17 0.1057178116    OMA17  RA5-A_18        A1        NbO-A5
## 187  Otu19  WGA54 0.1030440148    WGA54 WGA_YE5-B        B1       NbP-YE5
## 6   Others   WGA4 0.1013986014     WGA4 WGA_RA5-A        A1        NbO-A5
## 330   Otu5  WGA13 0.1001645413    WGA13 WGA_RA5-B        A1        NbO-A5
## 392   Otu7  WGA13 0.0993418346    WGA13 WGA_RA5-B        A1        NbO-A5
## 303   Otu4   WGA5 0.0964623612     WGA5 WGA_RA5-A        A1        NbO-A5
## 386   Otu7  WGA27 0.0917317976    WGA27 WGA_RD4-B        A2        NbO-D4
## 319   Otu5  OMA13 0.0878239408    OMA13  NB7-C_18        C1        NbQ-B7
## 363   Otu6  WGA34 0.0876182641    WGA34 WGA_NB7-A        C1        NbQ-B7
## 185  Otu19   WGA4 0.0841217606     WGA4 WGA_RA5-A        A1        NbO-A5
## 71   Otu11  OMA18 0.0839160839    OMA18  RA5-B_18        A1        NbO-A5
## 194  Otu19  WGA25 0.0830933772    WGA25 WGA_RD4-B        A2        NbO-D4
## 377   Otu7  OMA22 0.0818593172    OMA22  RD4-C_18        A2        NbO-D4
## 368   Otu6  OMA27 0.0812422871    OMA27  YE5-B_18        B1       NbP-YE5
## 352   Otu6  OMA28 0.0793911970    OMA28  YE5-C_18        B1       NbP-YE5
## 115  Otu12  WGA48 0.0781571370    WGA48 WGA_YE5-A        B1       NbP-YE5
## 191  Otu19   WGA5 0.0765117236     WGA5 WGA_RA5-A        A1        NbO-A5
## 365   Otu6  OMA13 0.0765117236    OMA13  NB7-C_18        C1        NbQ-B7
## 310   Otu5  OMA28 0.0758946935    OMA28  YE5-C_18        B1       NbP-YE5
## 221   Otu2  WGA13 0.0752776635    WGA13 WGA_RA5-B        A1        NbO-A5
## 373   Otu6  OMA20 0.0744549568    OMA20  RD4-A_18        A2        NbO-D4
## 280   Otu4  WGA45 0.0738379268    WGA45 WGA_YE5-A        B1       NbP-YE5
## 334   Otu5  WGA55 0.0736322501    WGA55 WGA_YE5-B        B1       NbP-YE5
## 229   Otu2  WGA46 0.0734265734    WGA46 WGA_YE5-A        B1       NbP-YE5
## 314   Otu5  OMA26 0.0709584533    OMA26  YE5-A_18        B1       NbP-YE5
## 146  Otu14  OMA11 0.0707527766    OMA11  NB7-A_18        C1        NbQ-B7
## 60   Otu10  WGA39 0.0699300699    WGA39 WGA_NB7-B        C1        NbQ-B7
## 324   Otu5  WGA27 0.0695187166    WGA27 WGA_RD4-B        A2        NbO-D4
## 344   Otu6  WGA20 0.0691073632    WGA20 WGA_RD4-A        A2        NbO-D4
## 336   Otu5  WGA53 0.0689016865    WGA53 WGA_YE5-B        B1       NbP-YE5
## 11  Others  WGA33 0.0682846565    WGA33 WGA_NB7-A        C1        NbQ-B7
## 378   Otu7  OMA21 0.0674619498    OMA21  RD4-B_18        A2        NbO-D4
## 341   Otu6  OMA11 0.0658165364    OMA11  NB7-A_18        C1        NbQ-B7
## 85   Otu11  WGA13 0.0656108597    WGA13 WGA_RA5-B        A1        NbO-A5
## 369   Otu6  OMA26 0.0656108597    OMA26  YE5-A_18        B1       NbP-YE5
## 103  Otu12  OMA27 0.0651995064    OMA27  YE5-B_18        B1       NbP-YE5
## 212   Otu2  WGA33 0.0651995064    WGA33 WGA_NB7-A        C1        NbQ-B7
## 338   Otu5  WGA48 0.0631427396    WGA48 WGA_YE5-A        B1       NbP-YE5
## 248  Otu24  WGA20 0.0625257096    WGA20 WGA_RD4-A        A2        NbO-D4
## 308   Otu5  OMA27 0.0625257096    OMA27  YE5-B_18        B1       NbP-YE5
## 37   Otu10  OMA26 0.0617030029    OMA26  YE5-A_18        B1       NbP-YE5
## 192  Otu19   WGA3 0.0573837927     WGA3 WGA_RA5-A        A1        NbO-A5
## 169  Otu14  WGA40 0.0557383793    WGA40 WGA_NB7-B        C1        NbQ-B7
## 371   Otu6  WGA39 0.0555327026    WGA39 WGA_NB7-B        C1        NbQ-B7
## 278   Otu4  WGA25 0.0547099959    WGA25 WGA_RD4-B        A2        NbO-D4
## 309   Otu5  OMA12 0.0547099959    OMA12  NB7-B_18        C1        NbQ-B7
## 298   Otu4  WGA12 0.0545043192    WGA12 WGA_RA5-B        A1        NbO-A5
## 302   Otu4  WGA40 0.0530645825    WGA40 WGA_NB7-B        C1        NbQ-B7
## 318   Otu5  OMA20 0.0526532291    OMA20  RD4-A_18        A2        NbO-D4
## 323   Otu5   WGA3 0.0514191691     WGA3 WGA_RA5-A        A1        NbO-A5
## 132  Otu12  WGA55 0.0501851090    WGA55 WGA_YE5-B        B1       NbP-YE5
## 23  Others  WGA34 0.0495680790    WGA34 WGA_NB7-A        C1        NbQ-B7
## 138  Otu14  WGA27 0.0495680790    WGA27 WGA_RD4-B        A2        NbO-D4
## 197  Otu19  OMA20 0.0485396956    OMA20  RD4-A_18        A2        NbO-D4
## 135  Otu12  WGA41 0.0481283422    WGA41 WGA_NB7-B        C1        NbQ-B7
## 307   Otu5  OMA11 0.0481283422    OMA11  NB7-A_18        C1        NbQ-B7
## 353   Otu6  OMA12 0.0479226656    OMA12  NB7-B_18        C1        NbQ-B7
## 53   Otu10  WGA46 0.0468942822    WGA46 WGA_YE5-A        B1       NbP-YE5
## 77   Otu11   WGA3 0.0462772522     WGA3 WGA_RA5-A        A1        NbO-A5
## 355   Otu6  WGA27 0.0450431921    WGA27 WGA_RD4-B        A2        NbO-D4
## 364   Otu6  WGA10 0.0444261621    WGA10 WGA_RA5-B        A1        NbO-A5
## 335   Otu5  WGA40 0.0431921020    WGA40 WGA_NB7-B        C1        NbQ-B7
## 96   Otu11  WGA10 0.0429864253    WGA10 WGA_RA5-B        A1        NbO-A5
## 31  Others  WGA19 0.0419580420    WGA19 WGA_RD4-A        A2        NbO-D4
## 38   Otu10  OMA12 0.0409296586    OMA12  NB7-B_18        C1        NbQ-B7
## 62   Otu10  WGA54 0.0407239819    WGA54 WGA_YE5-B        B1       NbP-YE5
## 36   Otu10  OMA11 0.0403126285    OMA11  NB7-A_18        C1        NbQ-B7
## 331   Otu5  OMA17 0.0399012752    OMA17  RA5-A_18        A1        NbO-A5
## 46   Otu10  OMA13 0.0394899218    OMA13  NB7-C_18        C1        NbQ-B7
## 345   Otu6  OMA21 0.0394899218    OMA21  RD4-B_18        A2        NbO-D4
## 19  Others  WGA18 0.0392842452    WGA18 WGA_RD4-A        A2        NbO-D4
## 315   Otu5  OMA19 0.0392842452    OMA19  RA5-C_18        A1        NbO-A5
## 282   Otu4  OMA19 0.0390785685    OMA19  RA5-C_18        A1        NbO-A5
## 81   Otu11  WGA20 0.0372274784    WGA20 WGA_RD4-A        A2        NbO-D4
## 40   Otu10  OMA27 0.0364047717    OMA27  YE5-B_18        B1       NbP-YE5
## 74   Otu11  OMA22 0.0361990950    OMA22  RD4-C_18        A2        NbO-D4
## 113  Otu12  OMA11 0.0359934183    OMA11  NB7-A_18        C1        NbQ-B7
## 275   Otu4  OMA27 0.0355820650    OMA27  YE5-B_18        B1       NbP-YE5
## 167  Otu14   WGA5 0.0347593583     WGA5 WGA_RA5-A        A1        NbO-A5
## 225   Otu2  WGA54 0.0347593583    WGA54 WGA_YE5-B        B1       NbP-YE5
## 48   Otu10  OMA28 0.0343480049    OMA28  YE5-C_18        B1       NbP-YE5
## 243  Otu24  OMA11 0.0339366516    OMA11  NB7-A_18        C1        NbQ-B7
## 159  Otu14  WGA19 0.0327025915    WGA19 WGA_RD4-A        A2        NbO-D4
## 252  Otu24  OMA28 0.0320855615    OMA28  YE5-C_18        B1       NbP-YE5
## 339   Otu5   WGA4 0.0320855615     WGA4 WGA_RA5-A        A1        NbO-A5
## 78   Otu11   WGA5 0.0312628548     WGA5 WGA_RA5-A        A1        NbO-A5
## 173  Otu19  OMA21 0.0312628548    OMA21  RD4-B_18        A2        NbO-D4
## 88   Otu11  OMA28 0.0296174414    OMA28  YE5-C_18        B1       NbP-YE5
## 140  Otu14  OMA28 0.0294117647    OMA28  YE5-C_18        B1       NbP-YE5
## 357   Otu6  OMA22 0.0294117647    OMA22  RD4-C_18        A2        NbO-D4
## 104  Otu12  OMA26 0.0290004114    OMA26  YE5-A_18        B1       NbP-YE5
## 41   Otu10  OMA20 0.0287947347    OMA20  RD4-A_18        A2        NbO-D4
## 154  Otu14  OMA12 0.0285890580    OMA12  NB7-B_18        C1        NbQ-B7
## 284   Otu4  OMA17 0.0285890580    OMA17  RA5-A_18        A1        NbO-A5
## 69   Otu11  OMA21 0.0281777046    OMA21  RD4-B_18        A2        NbO-D4
## 312   Otu5  OMA21 0.0279720280    OMA21  RD4-B_18        A2        NbO-D4
## 202  Otu19  WGA10 0.0277663513    WGA10 WGA_RA5-B        A1        NbO-A5
## 172  Otu19  OMA22 0.0275606746    OMA22  RD4-C_18        A2        NbO-D4
## 72   Otu11  OMA17 0.0271493213    OMA17  RA5-A_18        A1        NbO-A5
## 359   Otu6  WGA19 0.0263266146    WGA19 WGA_RD4-A        A2        NbO-D4
## 164  Otu14  WGA55 0.0248868778    WGA55 WGA_YE5-B        B1       NbP-YE5
## 316   Otu5  OMA18 0.0246812012    OMA18  RA5-B_18        A1        NbO-A5
## 325   Otu5  WGA25 0.0246812012    WGA25 WGA_RD4-B        A2        NbO-D4
## 283   Otu4  OMA18 0.0242698478    OMA18  RA5-B_18        A1        NbO-A5
## 190  Otu19  OMA28 0.0230357877    OMA28  YE5-C_18        B1       NbP-YE5
## 361   Otu6  OMA19 0.0230357877    OMA19  RA5-C_18        A1        NbO-A5
## 171  Otu19  OMA26 0.0220074044    OMA26  YE5-A_18        B1       NbP-YE5
## 52   Otu10  WGA48 0.0211846977    WGA48 WGA_YE5-A        B1       NbP-YE5
## 198  Otu19  OMA19 0.0203619910    OMA19  RA5-C_18        A1        NbO-A5
## 139  Otu14  WGA33 0.0199506376    WGA33 WGA_NB7-A        C1        NbQ-B7
## 178  Otu19  OMA27 0.0199506376    OMA27  YE5-B_18        B1       NbP-YE5
## 73   Otu11  OMA13 0.0191279309    OMA13  NB7-C_18        C1        NbQ-B7
## 125  Otu12  OMA12 0.0180995475    OMA12  NB7-B_18        C1        NbQ-B7
## 136  Otu12  WGA40 0.0180995475    WGA40 WGA_NB7-B        C1        NbQ-B7
## 261  Otu24  WGA13 0.0180995475    WGA13 WGA_RA5-B        A1        NbO-A5
## 273   Otu4  OMA26 0.0176881942    OMA26  YE5-A_18        B1       NbP-YE5
## 82   Otu11  OMA20 0.0174825175    OMA20  RD4-A_18        A2        NbO-D4
## 35   Otu10  OMA21 0.0172768408    OMA21  RD4-B_18        A2        NbO-D4
## 123  Otu12  WGA34 0.0170711641    WGA34 WGA_NB7-A        C1        NbQ-B7
## 79   Otu11  WGA27 0.0168654875    WGA27 WGA_RD4-B        A2        NbO-D4
## 76   Otu11  OMA27 0.0166598108    OMA27  YE5-B_18        B1       NbP-YE5
## 245  Otu24  OMA26 0.0166598108    OMA26  YE5-A_18        B1       NbP-YE5
## 372   Otu6  WGA12 0.0166598108    WGA12 WGA_RA5-B        A1        NbO-A5
## 83   Otu11  OMA19 0.0160427807    OMA19  RA5-C_18        A1        NbO-A5
## 264  Otu24  WGA10 0.0160427807    WGA10 WGA_RA5-B        A1        NbO-A5
## 137  Otu14   WGA3 0.0152200740     WGA3 WGA_RA5-A        A1        NbO-A5
## 161  Otu14  OMA19 0.0152200740    OMA19  RA5-C_18        A1        NbO-A5
## 204  Otu19  WGA12 0.0152200740    WGA12 WGA_RA5-B        A1        NbO-A5
## 347   Otu6   WGA4 0.0150143974     WGA4 WGA_RA5-A        A1        NbO-A5
## 63   Otu10  WGA53 0.0143973673    WGA53 WGA_YE5-B        B1       NbP-YE5
## 180  Otu19  WGA46 0.0139860140    WGA46 WGA_YE5-A        B1       NbP-YE5
## 51   Otu10  WGA27 0.0137803373    WGA27 WGA_RD4-B        A2        NbO-D4
## 265  Otu24  WGA34 0.0137803373    WGA34 WGA_NB7-A        C1        NbQ-B7
## 242  Otu24  OMA12 0.0133689840    OMA12  NB7-B_18        C1        NbQ-B7
## 260  Otu24  OMA18 0.0133689840    OMA18  RA5-B_18        A1        NbO-A5
## 327   Otu5  WGA45 0.0133689840    WGA45 WGA_YE5-A        B1       NbP-YE5
## 259  Otu24  OMA19 0.0125462773    OMA19  RA5-C_18        A1        NbO-A5
## 311   Otu5  OMA22 0.0119292472    OMA22  RD4-C_18        A2        NbO-D4
## 75   Otu11  OMA11 0.0117235705    OMA11  NB7-A_18        C1        NbQ-B7
## 241  Otu24  OMA13 0.0117235705    OMA13  NB7-C_18        C1        NbQ-B7
## 367   Otu6   WGA3 0.0117235705     WGA3 WGA_RA5-A        A1        NbO-A5
## 346   Otu6  WGA40 0.0115178939    WGA40 WGA_NB7-B        C1        NbQ-B7
## 244  Otu24  OMA27 0.0113122172    OMA27  YE5-B_18        B1       NbP-YE5
## 294   Otu4  OMA20 0.0113122172    OMA20  RD4-A_18        A2        NbO-D4
## 249  Otu24  OMA17 0.0111065405    OMA17  RA5-A_18        A1        NbO-A5
## 152  Otu14  WGA10 0.0106951872    WGA10 WGA_RA5-B        A1        NbO-A5
## 162  Otu14  OMA18 0.0104895105    OMA18  RA5-B_18        A1        NbO-A5
## 240  Otu24  OMA21 0.0104895105    OMA21  RD4-B_18        A2        NbO-D4
## 145  Otu14  OMA27 0.0094611271    OMA27  YE5-B_18        B1       NbP-YE5
## 263  Otu24  WGA12 0.0094611271    WGA12 WGA_RA5-B        A1        NbO-A5
## 349   Otu6  OMA18 0.0094611271    OMA18  RA5-B_18        A1        NbO-A5
## 366   Otu6   WGA5 0.0094611271     WGA5 WGA_RA5-A        A1        NbO-A5
## 374   Otu6  OMA17 0.0092554504    OMA17  RA5-A_18        A1        NbO-A5
## 342   Otu6  WGA46 0.0086384204    WGA46 WGA_YE5-A        B1       NbP-YE5
## 351   Otu6  WGA33 0.0086384204    WGA33 WGA_NB7-A        C1        NbQ-B7
## 175  Otu19  OMA13 0.0084327437    OMA13  NB7-C_18        C1        NbQ-B7
## 350   Otu6  WGA53 0.0082270671    WGA53 WGA_YE5-B        B1       NbP-YE5
## 43   Otu10  OMA18 0.0080213904    OMA18  RA5-B_18        A1        NbO-A5
## 200  Otu19  WGA18 0.0080213904    WGA18 WGA_RD4-A        A2        NbO-D4
## 250  Otu24  WGA54 0.0080213904    WGA54 WGA_YE5-B        B1       NbP-YE5
## 360   Otu6  WGA18 0.0080213904    WGA18 WGA_RD4-A        A2        NbO-D4
## 133  Otu12  WGA54 0.0078157137    WGA54 WGA_YE5-B        B1       NbP-YE5
## 124  Otu12  OMA13 0.0076100370    OMA13  NB7-C_18        C1        NbQ-B7
## 356   Otu6  WGA25 0.0076100370    WGA25 WGA_RD4-B        A2        NbO-D4
## 47   Otu10  WGA10 0.0074043603    WGA10 WGA_RA5-B        A1        NbO-A5
## 110  Otu12  OMA28 0.0074043603    OMA28  YE5-C_18        B1       NbP-YE5
## 116  Otu12  WGA46 0.0074043603    WGA46 WGA_YE5-A        B1       NbP-YE5
## 286   Otu4  OMA12 0.0074043603    OMA12  NB7-B_18        C1        NbQ-B7
## 112  Otu12   WGA3 0.0067873303     WGA3 WGA_RA5-A        A1        NbO-A5
## 274   Otu4  OMA11 0.0067873303    OMA11  NB7-A_18        C1        NbQ-B7
## 150  Otu14  OMA17 0.0065816536    OMA17  RA5-A_18        A1        NbO-A5
## 196  Otu19  WGA39 0.0065816536    WGA39 WGA_NB7-B        C1        NbQ-B7
## 55   Otu10  OMA19 0.0063759770    OMA19  RA5-C_18        A1        NbO-A5
## 45   Otu10  OMA17 0.0059646236    OMA17  RA5-A_18        A1        NbO-A5
## 156  Otu14  OMA26 0.0059646236    OMA26  YE5-A_18        B1       NbP-YE5
## 153  Otu14  OMA13 0.0057589469    OMA13  NB7-C_18        C1        NbQ-B7
## 160  Otu14  OMA20 0.0055532703    OMA20  RD4-A_18        A2        NbO-D4
## 254  Otu24   WGA3 0.0055532703     WGA3 WGA_RA5-A        A1        NbO-A5
## 287   Otu4  OMA28 0.0055532703    OMA28  YE5-C_18        B1       NbP-YE5
## 317   Otu5  WGA20 0.0055532703    WGA20 WGA_RD4-A        A2        NbO-D4
## 87   Otu11  OMA12 0.0049362402    OMA12  NB7-B_18        C1        NbQ-B7
## 234   Otu2  WGA41 0.0049362402    WGA41 WGA_NB7-B        C1        NbQ-B7
## 42   Otu10  OMA22 0.0047305636    OMA22  RD4-C_18        A2        NbO-D4
## 177  Otu19  OMA11 0.0047305636    OMA11  NB7-A_18        C1        NbQ-B7
## 182  Otu19  WGA20 0.0047305636    WGA20 WGA_RD4-A        A2        NbO-D4
## 375   Otu7  OMA11 0.0045248869    OMA11  NB7-A_18        C1        NbQ-B7
## 176  Otu19  OMA12 0.0043192102    OMA12  NB7-B_18        C1        NbQ-B7
## 404   Otu7  WGA18 0.0043192102    WGA18 WGA_RD4-A        A2        NbO-D4
## 102  Otu11  WGA55 0.0041135335    WGA55 WGA_YE5-B        B1       NbP-YE5
## 70   Otu11  OMA26 0.0037021802    OMA26  YE5-A_18        B1       NbP-YE5
## 158  Otu14  OMA21 0.0037021802    OMA21  RD4-B_18        A2        NbO-D4
## 383   Otu7  OMA27 0.0037021802    OMA27  YE5-B_18        B1       NbP-YE5
## 328   Otu5  WGA19 0.0034965035    WGA19 WGA_RD4-A        A2        NbO-D4
## 271  Otu24   WGA4 0.0032908268     WGA4 WGA_RA5-A        A1        NbO-A5
## 320   Otu5  WGA10 0.0032908268    WGA10 WGA_RA5-B        A1        NbO-A5
## 385   Otu7  OMA26 0.0032908268    OMA26  YE5-A_18        B1       NbP-YE5
## 239  Otu24  OMA22 0.0028794735    OMA22  RD4-C_18        A2        NbO-D4
## 108  Otu12  OMA20 0.0024681201    OMA20  RD4-A_18        A2        NbO-D4
## 258  Otu24  OMA20 0.0024681201    OMA20  RD4-A_18        A2        NbO-D4
## 381   Otu7  OMA12 0.0024681201    OMA12  NB7-B_18        C1        NbQ-B7
## 157  Otu14  OMA22 0.0022624434    OMA22  RD4-C_18        A2        NbO-D4
## 119  Otu12  OMA18 0.0018510901    OMA18  RA5-B_18        A1        NbO-A5
## 120  Otu12  OMA17 0.0016454134    OMA17  RA5-A_18        A1        NbO-A5
## 285   Otu4  OMA13 0.0014397367    OMA13  NB7-C_18        C1        NbQ-B7
## 380   Otu7  OMA13 0.0014397367    OMA13  NB7-C_18        C1        NbQ-B7
## 67   Otu10  WGA19 0.0012340601    WGA19 WGA_RD4-A        A2        NbO-D4
## 80   Otu11  WGA25 0.0012340601    WGA25 WGA_RD4-B        A2        NbO-D4
## 109  Otu12  OMA22 0.0012340601    OMA22  RD4-C_18        A2        NbO-D4
## 118  Otu12  OMA19 0.0012340601    OMA19  RA5-C_18        A1        NbO-A5
## 170  Otu14   WGA4 0.0012340601     WGA4 WGA_RA5-A        A1        NbO-A5
## 193  Otu19  WGA27 0.0012340601    WGA27 WGA_RD4-B        A2        NbO-D4
## 382   Otu7  OMA28 0.0012340601    OMA28  YE5-C_18        B1       NbP-YE5
## 255  Otu24   WGA5 0.0010283834     WGA5 WGA_RA5-A        A1        NbO-A5
## 105  Otu12  OMA21 0.0008227067    OMA21  RD4-B_18        A2        NbO-D4
## 389   Otu7  OMA20 0.0008227067    OMA20  RD4-A_18        A2        NbO-D4
## 57   Otu10  WGA13 0.0006170300    WGA13 WGA_RA5-B        A1        NbO-A5
## 304   Otu4  WGA39 0.0006170300    WGA39 WGA_NB7-B        C1        NbQ-B7
## 58   Otu10   WGA4 0.0004113534     WGA4 WGA_RA5-A        A1        NbO-A5
## 266  Otu24  WGA40 0.0004113534    WGA40 WGA_NB7-B        C1        NbQ-B7
## 272  Otu24  WGA39 0.0004113534    WGA39 WGA_NB7-B        C1        NbQ-B7
## 329   Otu5  WGA18 0.0004113534    WGA18 WGA_RD4-A        A2        NbO-D4
## 49   Otu10  WGA33 0.0002056767    WGA33 WGA_NB7-A        C1        NbQ-B7
## 64   Otu10   WGA5 0.0002056767     WGA5 WGA_RA5-A        A1        NbO-A5
## 65   Otu10  WGA41 0.0002056767    WGA41 WGA_NB7-B        C1        NbQ-B7
## 126  Otu12   WGA5 0.0002056767     WGA5 WGA_RA5-A        A1        NbO-A5
## 267  Otu24  WGA45 0.0002056767    WGA45 WGA_YE5-A        B1       NbP-YE5
## 348   Otu6  WGA13 0.0002056767    WGA13 WGA_RA5-B        A1        NbO-A5
## 39   Otu10  WGA25 0.0000000000    WGA25 WGA_RD4-B        A2        NbO-D4
## 44   Otu10  WGA20 0.0000000000    WGA20 WGA_RD4-A        A2        NbO-D4
## 50   Otu10   WGA3 0.0000000000     WGA3 WGA_RA5-A        A1        NbO-A5
## 54   Otu10  WGA45 0.0000000000    WGA45 WGA_YE5-A        B1       NbP-YE5
## 56   Otu10  WGA18 0.0000000000    WGA18 WGA_RD4-A        A2        NbO-D4
## 59   Otu10  WGA12 0.0000000000    WGA12 WGA_RA5-B        A1        NbO-A5
## 61   Otu10  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 66   Otu10  WGA40 0.0000000000    WGA40 WGA_NB7-B        C1        NbQ-B7
## 68   Otu10  WGA55 0.0000000000    WGA55 WGA_YE5-B        B1       NbP-YE5
## 89   Otu11  WGA33 0.0000000000    WGA33 WGA_NB7-A        C1        NbQ-B7
## 90   Otu11  WGA53 0.0000000000    WGA53 WGA_YE5-B        B1       NbP-YE5
## 91   Otu11  WGA45 0.0000000000    WGA45 WGA_YE5-A        B1       NbP-YE5
## 92   Otu11  WGA19 0.0000000000    WGA19 WGA_RD4-A        A2        NbO-D4
## 93   Otu11  WGA48 0.0000000000    WGA48 WGA_YE5-A        B1       NbP-YE5
## 94   Otu11  WGA46 0.0000000000    WGA46 WGA_YE5-A        B1       NbP-YE5
## 97   Otu11  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
## 100  Otu11  WGA54 0.0000000000    WGA54 WGA_YE5-B        B1       NbP-YE5
## 101  Otu11  WGA39 0.0000000000    WGA39 WGA_NB7-B        C1        NbQ-B7
## 106  Otu12  WGA25 0.0000000000    WGA25 WGA_RD4-B        A2        NbO-D4
## 107  Otu12  WGA20 0.0000000000    WGA20 WGA_RD4-A        A2        NbO-D4
## 114  Otu12  WGA27 0.0000000000    WGA27 WGA_RD4-B        A2        NbO-D4
## 121  Otu12  WGA12 0.0000000000    WGA12 WGA_RA5-B        A1        NbO-A5
## 122  Otu12  WGA10 0.0000000000    WGA10 WGA_RA5-B        A1        NbO-A5
## 127  Otu12  WGA19 0.0000000000    WGA19 WGA_RD4-A        A2        NbO-D4
## 128  Otu12  WGA18 0.0000000000    WGA18 WGA_RD4-A        A2        NbO-D4
## 129  Otu12  WGA13 0.0000000000    WGA13 WGA_RA5-B        A1        NbO-A5
## 130  Otu12   WGA4 0.0000000000     WGA4 WGA_RA5-A        A1        NbO-A5
## 142  Otu14  WGA20 0.0000000000    WGA20 WGA_RD4-A        A2        NbO-D4
## 143  Otu14  WGA46 0.0000000000    WGA46 WGA_YE5-A        B1       NbP-YE5
## 144  Otu14  WGA45 0.0000000000    WGA45 WGA_YE5-A        B1       NbP-YE5
## 147  Otu14  WGA18 0.0000000000    WGA18 WGA_RD4-A        A2        NbO-D4
## 149  Otu14  WGA12 0.0000000000    WGA12 WGA_RA5-B        A1        NbO-A5
## 151  Otu14  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 155  Otu14  WGA48 0.0000000000    WGA48 WGA_YE5-A        B1       NbP-YE5
## 163  Otu14  WGA39 0.0000000000    WGA39 WGA_NB7-B        C1        NbQ-B7
## 165  Otu14  WGA54 0.0000000000    WGA54 WGA_YE5-B        B1       NbP-YE5
## 166  Otu14  WGA53 0.0000000000    WGA53 WGA_YE5-B        B1       NbP-YE5
## 179  Otu19  WGA48 0.0000000000    WGA48 WGA_YE5-A        B1       NbP-YE5
## 183  Otu19  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
## 184  Otu19  WGA40 0.0000000000    WGA40 WGA_NB7-B        C1        NbQ-B7
## 188  Otu19  WGA53 0.0000000000    WGA53 WGA_YE5-B        B1       NbP-YE5
## 189  Otu19  WGA33 0.0000000000    WGA33 WGA_NB7-A        C1        NbQ-B7
## 195  Otu19  WGA55 0.0000000000    WGA55 WGA_YE5-B        B1       NbP-YE5
## 201  Otu19  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 246  Otu24  WGA27 0.0000000000    WGA27 WGA_RD4-B        A2        NbO-D4
## 247  Otu24  WGA25 0.0000000000    WGA25 WGA_RD4-B        A2        NbO-D4
## 251  Otu24  WGA53 0.0000000000    WGA53 WGA_YE5-B        B1       NbP-YE5
## 253  Otu24  WGA33 0.0000000000    WGA33 WGA_NB7-A        C1        NbQ-B7
## 257  Otu24  WGA46 0.0000000000    WGA46 WGA_YE5-A        B1       NbP-YE5
## 262  Otu24  WGA55 0.0000000000    WGA55 WGA_YE5-B        B1       NbP-YE5
## 268  Otu24  WGA19 0.0000000000    WGA19 WGA_RD4-A        A2        NbO-D4
## 269  Otu24  WGA18 0.0000000000    WGA18 WGA_RD4-A        A2        NbO-D4
## 270  Otu24  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
## 291   Otu4  WGA48 0.0000000000    WGA48 WGA_YE5-A        B1       NbP-YE5
## 300   Otu4  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 301   Otu4  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
## 321   Otu5  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 333   Otu5  WGA39 0.0000000000    WGA39 WGA_NB7-B        C1        NbQ-B7
## 340   Otu5  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
## 376   Otu7  WGA20 0.0000000000    WGA20 WGA_RD4-A        A2        NbO-D4
## 388   Otu7  WGA46 0.0000000000    WGA46 WGA_YE5-A        B1       NbP-YE5
## 395   Otu7  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 396   Otu7  WGA33 0.0000000000    WGA33 WGA_NB7-A        C1        NbQ-B7
## 397   Otu7  WGA53 0.0000000000    WGA53 WGA_YE5-B        B1       NbP-YE5
## 399   Otu7  WGA19 0.0000000000    WGA19 WGA_RD4-A        A2        NbO-D4
## 400   Otu7  WGA48 0.0000000000    WGA48 WGA_YE5-A        B1       NbP-YE5
## 401   Otu7  WGA40 0.0000000000    WGA40 WGA_NB7-B        C1        NbQ-B7
## 402   Otu7  WGA45 0.0000000000    WGA45 WGA_YE5-A        B1       NbP-YE5
## 403   Otu7  WGA39 0.0000000000    WGA39 WGA_NB7-B        C1        NbQ-B7
## 405   Otu7  WGA54 0.0000000000    WGA54 WGA_YE5-B        B1       NbP-YE5
## 407   Otu7  WGA55 0.0000000000    WGA55 WGA_YE5-B        B1       NbP-YE5
## 408   Otu7  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
##     Treatment Treatment3 Treatment4      Pop Pop2 Pop3 Sample_Site
## 296       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 295       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 281       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 111       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 18        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 224       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 231       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 222       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 7         NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 279       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 214       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 358       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 390       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 205       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 207       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 289       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 362       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 206       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 406       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 299       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 12        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 290       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 379       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 211       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 34        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 141       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 220       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 210       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 99        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 326       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 384       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 292       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 5         NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 13        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 256       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 305       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 181       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 277       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 28        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 1         RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 2         YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 30        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 32        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 168       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 21        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 235       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 4         YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 10        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 86        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 25        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 20        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 387       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 208       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 218       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 236       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 14        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 26        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 219       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 398       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 9         RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 288       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 24        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 276       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 228       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 223       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 22        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 3         RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 199       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 332       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 17        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 16        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 216       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 337       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 230       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 33        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 226       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 203       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 354       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 84        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 117       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 215       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 27        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 293       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 98        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 233       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 8         RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 131       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 306       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 186       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 29        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 213       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 134       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 95        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 238       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 322       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 393       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 391       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 343       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 227       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 297       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 15        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 148       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 237       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 370       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 217       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 394       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 209       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 313       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 232       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 174       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 187       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 6         RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 330       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 392       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 303       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 386       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 319       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 363       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 185       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 71        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 194       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 377       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 368       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 352       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 115       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 191       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 365       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 310       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 221       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 373       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 280       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 334       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 229       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 314       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 146       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 60        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 324       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 344       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 336       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 11        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 378       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 341       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 85        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 369       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 103       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 212       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 338       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 248       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 308       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 37        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 192       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 169       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 371       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 278       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 309       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 298       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 302       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 318       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 323       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 132       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 23        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 138       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 197       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 135       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 307       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 353       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 53        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 77        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 355       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 364       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 335       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 96        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 31        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 38        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 62        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 36        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 331       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 46        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 345       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 19        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 315       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 282       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 81        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 40        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 74        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 113       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 275       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 167       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 225       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 48        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 243       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 159       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 252       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 339       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 78        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 173       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 88        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 140       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 357       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 104       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 41        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 154       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 284       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 69        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 312       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 202       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 172       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 72        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 359       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 164       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 316       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 325       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 283       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 190       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 361       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 171       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 52        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 198       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 139       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 178       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 73        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 125       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 136       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 261       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 273       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 82        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 35        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 123       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 79        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 76        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 245       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 372       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 83        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 264       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 137       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 161       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 204       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 347       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 63        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 180       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 51        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 265       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 242       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 260       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 327       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 259       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 311       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 75        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 241       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 367       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 346       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 244       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 294       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 249       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 152       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 162       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 240       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 145       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 263       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 349       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 366       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 374       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 342       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 351       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 175       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 350       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 43        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 200       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 250       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 360       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 133       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 124       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 356       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 47        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 110       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 116       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 286       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 112       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 274       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 150       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 196       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 55        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 45        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 156       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 153       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 160       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 254       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 287       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 317       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 87        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 234       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 42        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 177       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 182       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 375       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 176       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 404       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 102       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 70        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 158       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 383       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 328       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 271       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 320       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 385       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 239       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 108       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 258       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 381       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 157       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 119       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 120       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 285       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 380       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 67        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 80        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 109       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 118       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 170       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 193       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 382       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 255       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 105       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 389       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 57        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 304       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 58        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 266       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 272       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 329       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 49        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 64        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 65        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 126       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 267       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 348       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 39        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 44        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 50        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 54        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 56        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 59        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 61        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 66        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 68        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 89        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 90        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 91        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 92        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 93        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 94        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 97        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 100       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 101       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 106       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 107       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 114       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 121       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 122       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 127       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 128       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 129       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 130       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 142       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 143       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 144       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 147       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 149       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 151       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 155       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 163       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 165       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 166       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 179       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 183       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 184       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 188       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 189       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 195       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 201       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 246       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 247       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 251       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 253       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 257       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 262       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 268       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 269       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 270       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 291       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 300       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 301       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 321       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 333       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 340       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 376       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 388       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 395       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 396       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 397       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 399       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 400       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 401       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 402       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 403       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 405       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 407       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 408       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
##      treatment_26           treatment_26.1 extract         definition_26 temp
## 296  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 295  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 281     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 111  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 18   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 224  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 231     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 222     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 7    NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 279  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 214     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 358  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 390     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 205    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 207    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 289 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 362 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 206     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 406  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 299  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 12     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 290  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 379     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 211    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 34  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 141  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 220     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 210     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 99   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 326 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 384  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 292 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 5       NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 13      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 256 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 305 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 181 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 277     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 28     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 1    NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 2   NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 30   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 32      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 168  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 21      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 235  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 4      NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 10  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 86   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 25      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 20   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 387  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 208     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 218  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 236  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 14  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 26  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 219     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 398  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 9       NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 288  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 24   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 276  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 228 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 223  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 22  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 3    NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 199  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 332  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 17      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 16   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 216  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 337  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 230  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 33      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 226 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 203  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 354 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 84   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 117 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 215  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 27   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 293  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 98   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 233  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 8    NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 131  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 306  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 186     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 29      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 213  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 134 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 95   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 238  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 322 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 393  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 391     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 343 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 227  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 297 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 15   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 148  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 237 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 370 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 217 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 394  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 209     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 313  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 232  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 174     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 187 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 6    NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 330  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 392  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 303  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 386  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 319     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 363  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 185  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 71      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 194  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 377     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 368    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 352    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 115 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 191  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 365     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 310    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 221  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 373     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 280 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 334 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 229 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 314    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 146     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 60   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 324  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 344  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 336 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 11   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 378     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 341     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 85   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 369    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 103    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 212  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 338 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 248  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 308    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 37     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 192  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 169  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 371  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 278  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 309     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 298  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 302  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 318     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 323  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 132 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 23   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 138  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 197     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 135  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 307     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 353     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 53  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 77   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 355  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 364  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 335  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 96   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 31   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 38      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 62  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 36      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 331     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 46      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 345     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 19   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 315     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 282     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 81   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 40     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 74      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 113     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 275    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 167  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 225 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 48     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 243     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 159  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 252    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 339  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 78   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 173     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 88     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 140    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 357     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 104    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 41      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 154     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 284     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 69      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 312     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 202  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 172     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 72      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 359  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 164 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 316     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 325  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 283     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 190    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 361     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 171    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 52  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 198     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 139  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 178    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 73      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 125     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 136  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 261  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 273    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 82      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 35      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 123  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 79   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 76     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 245    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 372  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 83      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 264  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 137  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 161     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 204  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 347  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 63  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 180 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 51   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 265  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 242     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 260     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 327 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 259     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 311     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 75      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 241     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 367  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 346  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 244    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 294     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 249     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 152  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 162     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 240     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 145    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 263  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 349     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 366  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 374     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 342 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 351  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 175     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 350 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 43      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 200  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 250 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 360  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 133 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 124     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 356  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 47   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 110    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 116 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 286     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 112  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 274     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 150     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 196  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 55      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 45      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 156    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 153     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 160     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 254  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 287    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 317  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 87      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 234  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 42      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 177     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 182  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 375     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 176     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 404  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 102 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 70     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 158     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 383    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 328  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 271  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 320  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 385    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 239     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 108     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 258     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 381     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 157     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 119     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 120     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 285     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 380     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 67   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 80   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 109     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 118     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 170  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 193  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 382    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 255  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 105     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 389     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 57   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 304  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 58   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 266  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 272  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 329  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 49   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 64   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 65   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 126  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 267 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 348  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 39   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 44   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 50   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 54  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 56   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 59   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 61   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 66   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 68  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 89   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 90  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 91  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 92   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 93  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 94  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 97   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 100 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 101  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 106  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 107  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 114  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 121  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 122  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 127  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 128  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 129  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 130  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 142  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 143 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 144 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 147  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 149  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 151  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 155 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 163  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 165 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 166 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 179 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 183  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 184  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 188 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 189  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 195 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 201  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 246  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 247  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 251 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 253  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 257 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 262 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 268  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 269  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 270  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 291 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 300  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 301  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 321  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 333  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 340  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 376  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 388 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 395  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 396  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 397 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 399  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 400 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 401  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 402 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 403  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 405 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 407 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 408  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
##             col    col2  col_26 pch culture definition keep_2026 diatom_control
## 296        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 295        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 281        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 111 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 18  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 224 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 231        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 222 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 7   forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 279        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 214 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 358 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 390        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 205   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 207   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 289   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 362   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 206 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 406  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 299  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 12    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 290        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 379        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 211   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 34    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 141        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 220        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 210        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 99  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 326   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 384  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 292   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 5   forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 13  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 256   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 305   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 181   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 277        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 28    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 1    steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 2     firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 30  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 32         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 168 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 21         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 235  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 4     firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 10    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 86   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 25  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 20   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 387        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 208        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 218        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 236 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 14    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 26    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 219        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 398  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 9          navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 288 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 24   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 276  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 228   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 223  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 22    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 3    steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 199        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 332  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 17         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 16         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 216        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 337  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 230 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 33         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 226   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 203  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 354   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 84         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 117   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 215        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 27         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 293  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 98  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 233        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 8    steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 131 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 306  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 186        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 29         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 213  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 134   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 95   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 238  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 322   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 393  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 391        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 343   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 227  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 297   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 15         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 148  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 237   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 370   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 217   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 394  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 209        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 313 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 232        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 174        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 187   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 6    steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 330  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 392  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 303  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 386        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 319 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 363 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 185  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 71         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 194        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 377        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 368   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 352   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 115   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 191  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 365 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 310   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 221  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 373        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 280   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 334   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 229   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 314   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 146 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 60  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 324        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 344        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 336   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 11  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 378        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 341 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 85   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 369   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 103   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 212 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 338   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 248        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 308   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 37    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 192  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 169 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 371 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 278        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 309 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 298  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 302 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 318        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 323  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 132   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 23  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 138        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 197        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 135 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 307 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 353 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 53    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 77   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 355        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 364  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 335 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 96   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 31         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 38  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 62    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 36  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 331        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 46  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 345        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 19         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 315        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 282        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 81         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 40    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 74         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 113 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 275   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 167  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 225   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 48    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 243 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 159        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 252   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 339  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 78   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 173        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 88    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 140   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 357        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 104   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 41         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 154 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 284        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 69         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 312        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 202  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 172        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 72         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 359        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 164   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 316        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 325        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 283        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 190   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 361        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 171   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 52    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 198        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 139 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 178   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 73  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 125 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 136 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 261  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 273   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 82         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 35         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 123 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 79         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 76    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 245   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 372  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 83         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 264  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 137  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 161        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 204  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 347  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 63    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 180   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 51         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 265 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 242 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 260        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 327   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 259        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 311        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 75  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 241 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 367  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 346 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 244   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 294        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 249        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 152  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 162        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 240        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 145   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 263  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 349        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 366  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 374        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 342   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 351 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 175 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 350   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 43         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 200        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 250   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 360        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 133   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 124 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 356        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 47   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 110   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 116   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 286 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 112  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 274 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 150        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 196 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 55         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 45         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 156   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 153 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 160        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 254  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 287   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 317        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 87  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 234 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 42         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 177 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 182        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 375 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 176 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 404        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 102   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 70    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 158        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 383   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 328        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 271  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 320  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 385   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 239        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 108        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 258        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 381 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 157        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 119        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 120        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 285 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 380 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 67         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 80         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 109        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 118        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 170  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 193        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 382   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 255  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 105        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 389        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 57   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 304 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 58   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 266 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 272 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 329        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 49  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 64   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 65  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 126  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 267   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 348  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 39         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 44         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 50   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 54    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 56         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 59   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 61  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 66  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 68    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 89  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 90    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 91    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 92         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 93    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 94    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 97  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 100   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 101 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 106        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 107        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 114        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 121  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 122  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 127        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 128        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 129  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 130  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 142        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 143   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 144   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 147        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 149  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 151 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 155   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 163 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 165   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 166   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 179   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 183 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 184 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 188   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 189 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 195   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 201 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 246        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 247        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 251   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 253 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 257   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 262   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 268        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 269        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 270 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 291   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 300 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 301 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 321 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 333 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 340 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 376        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 388   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 395 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 396 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 397   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 399        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 400   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 401 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 402   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 403 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 405   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 407   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 408 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
##     innoc  Kingdom          phyla               class            family
## 296  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 295  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 281  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 111  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 18   <NA>     <NA>           <NA>                <NA>              <NA>
## 224  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 231  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 222  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 7    <NA>     <NA>           <NA>                <NA>              <NA>
## 279  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 214  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 358  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 390  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 205  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 207  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 289  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 362  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 206  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 406  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 299  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 12   <NA>     <NA>           <NA>                <NA>              <NA>
## 290  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 379  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 211  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 34   <NA>     <NA>           <NA>                <NA>              <NA>
## 141  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 220  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 210  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 99   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 326  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 384  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 292  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 5    <NA>     <NA>           <NA>                <NA>              <NA>
## 13   <NA>     <NA>           <NA>                <NA>              <NA>
## 256  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 305  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 181  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 277  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 28   <NA>     <NA>           <NA>                <NA>              <NA>
## 1    <NA>     <NA>           <NA>                <NA>              <NA>
## 2    <NA>     <NA>           <NA>                <NA>              <NA>
## 30   <NA>     <NA>           <NA>                <NA>              <NA>
## 32   <NA>     <NA>           <NA>                <NA>              <NA>
## 168  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 21   <NA>     <NA>           <NA>                <NA>              <NA>
## 235  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 4    <NA>     <NA>           <NA>                <NA>              <NA>
## 10   <NA>     <NA>           <NA>                <NA>              <NA>
## 86   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 25   <NA>     <NA>           <NA>                <NA>              <NA>
## 20   <NA>     <NA>           <NA>                <NA>              <NA>
## 387  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 208  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 218  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 236  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 14   <NA>     <NA>           <NA>                <NA>              <NA>
## 26   <NA>     <NA>           <NA>                <NA>              <NA>
## 219  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 398  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 9    <NA>     <NA>           <NA>                <NA>              <NA>
## 288  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 24   <NA>     <NA>           <NA>                <NA>              <NA>
## 276  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 228  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 223  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 22   <NA>     <NA>           <NA>                <NA>              <NA>
## 3    <NA>     <NA>           <NA>                <NA>              <NA>
## 199  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 332  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 17   <NA>     <NA>           <NA>                <NA>              <NA>
## 16   <NA>     <NA>           <NA>                <NA>              <NA>
## 216  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 337  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 230  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 33   <NA>     <NA>           <NA>                <NA>              <NA>
## 226  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 203  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 354  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 84   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 117  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 215  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 27   <NA>     <NA>           <NA>                <NA>              <NA>
## 293  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 98   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 233  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 8    <NA>     <NA>           <NA>                <NA>              <NA>
## 131  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 306  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 186  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 29   <NA>     <NA>           <NA>                <NA>              <NA>
## 213  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 134  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 95   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 238  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 322  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 393  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 391  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 343  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 227  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 297  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 15   <NA>     <NA>           <NA>                <NA>              <NA>
## 148  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 237  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 370  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 217  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 394  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 209  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 313  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 232  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 174  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 187  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 6    <NA>     <NA>           <NA>                <NA>              <NA>
## 330  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 392  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 303  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 386  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 319  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 363  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 185  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 71   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 194  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 377  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 368  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 352  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 115  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 191  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 365  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 310  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 221  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 373  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 280  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 334  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 229  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 314  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 146  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 60   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 324  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 344  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 336  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 11   <NA>     <NA>           <NA>                <NA>              <NA>
## 378  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 341  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 85   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 369  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 103  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 212  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 338  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 248  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 308  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 37   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 192  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 169  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 371  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 278  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 309  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 298  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 302  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 318  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 323  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 132  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 23   <NA>     <NA>           <NA>                <NA>              <NA>
## 138  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 197  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 135  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 307  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 353  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 53   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 77   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 355  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 364  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 335  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 96   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 31   <NA>     <NA>           <NA>                <NA>              <NA>
## 38   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 62   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 36   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 331  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 46   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 345  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 19   <NA>     <NA>           <NA>                <NA>              <NA>
## 315  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 282  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 81   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 40   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 74   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 113  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 275  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 167  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 225  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 48   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 243  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 159  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 252  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 339  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 78   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 173  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 88   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 140  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 357  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 104  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 41   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 154  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 284  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 69   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 312  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 202  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 172  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 72   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 359  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 164  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 316  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 325  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 283  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 190  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 361  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 171  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 52   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 198  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 139  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 178  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 73   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 125  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 136  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 261  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 273  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 82   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 35   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 123  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 79   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 76   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 245  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 372  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 83   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 264  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 137  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 161  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 204  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 347  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 63   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 180  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 51   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 265  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 242  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 260  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 327  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 259  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 311  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 75   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 241  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 367  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 346  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 244  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 294  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 249  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 152  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 162  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 240  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 145  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 263  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 349  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 366  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 374  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 342  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 351  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 175  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 350  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 43   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 200  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 250  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 360  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 133  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 124  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 356  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 47   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 110  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 116  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 286  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 112  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 274  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 150  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 196  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 55   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 45   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 156  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 153  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 160  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 254  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 287  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 317  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 87   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 234  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 42   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 177  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 182  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 375  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 176  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 404  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 102  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 70   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 158  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 383  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 328  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 271  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 320  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 385  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 239  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 108  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 258  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 381  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 157  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 119  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 120  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 285  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 380  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 67   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 80   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 109  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 118  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 170  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 193  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 382  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 255  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 105  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 389  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 57   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 304  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 58   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 266  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 272  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 329  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 49   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 64   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 65   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 126  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 267  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 348  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 39   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 44   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 50   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 54   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 56   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 59   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 61   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 66   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 68   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 89   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 90   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 91   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 92   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 93   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 94   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 97   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 100  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 101  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 106  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 107  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 114  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 121  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 122  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 127  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 128  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 129  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 130  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 142  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 143  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 144  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 147  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 149  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 151  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 155  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 163  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 165  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 166  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 179  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 183  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 184  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 188  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 189  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 195  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 201  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 246  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 247  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 251  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 253  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 257  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 262  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 268  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 269  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 270  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 291  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 300  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 301  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 321  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 333  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 340  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 376  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 388  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 395  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 396  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 397  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 399  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 400  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 401  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 402  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 403  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 405  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 407  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 408  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
##                  genus                 species
## 296 Oceanospirillaceae             Marinomonas
## 295 Oceanospirillaceae             Marinomonas
## 281 Oceanospirillaceae             Marinomonas
## 111 Oceanospirillaceae          Neptuniibacter
## 18                <NA>                    <NA>
## 224  Flavobacteriaceae           Tenacibaculum
## 231  Flavobacteriaceae           Tenacibaculum
## 222  Flavobacteriaceae           Tenacibaculum
## 7                 <NA>                    <NA>
## 279 Oceanospirillaceae             Marinomonas
## 214  Flavobacteriaceae           Tenacibaculum
## 358              SAR11 Candidatus Pelagibacter
## 390       Vibrionaceae                  Vibrio
## 205  Flavobacteriaceae           Tenacibaculum
## 207  Flavobacteriaceae           Tenacibaculum
## 289 Oceanospirillaceae             Marinomonas
## 362              SAR11 Candidatus Pelagibacter
## 206  Flavobacteriaceae           Tenacibaculum
## 406       Vibrionaceae                  Vibrio
## 299 Oceanospirillaceae             Marinomonas
## 12                <NA>                    <NA>
## 290 Oceanospirillaceae             Marinomonas
## 379       Vibrionaceae                  Vibrio
## 211  Flavobacteriaceae           Tenacibaculum
## 34                <NA>                    <NA>
## 141 Oceanospirillaceae         Marinobacterium
## 220  Flavobacteriaceae           Tenacibaculum
## 210  Flavobacteriaceae           Tenacibaculum
## 99      Cryomorphaceae             Owenweeksia
## 326   Alteromonadaceae         Aestuariibacter
## 384       Vibrionaceae                  Vibrio
## 292 Oceanospirillaceae             Marinomonas
## 5                 <NA>                    <NA>
## 13                <NA>                    <NA>
## 256  Flavobacteriaceae              Aquibacter
## 305 Oceanospirillaceae             Marinomonas
## 181     Cryomorphaceae                Wandonia
## 277 Oceanospirillaceae             Marinomonas
## 28                <NA>                    <NA>
## 1                 <NA>                    <NA>
## 2                 <NA>                    <NA>
## 30                <NA>                    <NA>
## 32                <NA>                    <NA>
## 168 Oceanospirillaceae         Marinobacterium
## 21                <NA>                    <NA>
## 235  Flavobacteriaceae           Tenacibaculum
## 4                 <NA>                    <NA>
## 10                <NA>                    <NA>
## 86      Cryomorphaceae             Owenweeksia
## 25                <NA>                    <NA>
## 20                <NA>                    <NA>
## 387       Vibrionaceae                  Vibrio
## 208  Flavobacteriaceae           Tenacibaculum
## 218  Flavobacteriaceae           Tenacibaculum
## 236  Flavobacteriaceae           Tenacibaculum
## 14                <NA>                    <NA>
## 26                <NA>                    <NA>
## 219  Flavobacteriaceae           Tenacibaculum
## 398       Vibrionaceae                  Vibrio
## 9                 <NA>                    <NA>
## 288 Oceanospirillaceae             Marinomonas
## 24                <NA>                    <NA>
## 276 Oceanospirillaceae             Marinomonas
## 228  Flavobacteriaceae           Tenacibaculum
## 223  Flavobacteriaceae           Tenacibaculum
## 22                <NA>                    <NA>
## 3                 <NA>                    <NA>
## 199     Cryomorphaceae                Wandonia
## 332   Alteromonadaceae         Aestuariibacter
## 17                <NA>                    <NA>
## 16                <NA>                    <NA>
## 216  Flavobacteriaceae           Tenacibaculum
## 337   Alteromonadaceae         Aestuariibacter
## 230  Flavobacteriaceae           Tenacibaculum
## 33                <NA>                    <NA>
## 226  Flavobacteriaceae           Tenacibaculum
## 203     Cryomorphaceae                Wandonia
## 354              SAR11 Candidatus Pelagibacter
## 84      Cryomorphaceae             Owenweeksia
## 117 Oceanospirillaceae          Neptuniibacter
## 215  Flavobacteriaceae           Tenacibaculum
## 27                <NA>                    <NA>
## 293 Oceanospirillaceae             Marinomonas
## 98      Cryomorphaceae             Owenweeksia
## 233  Flavobacteriaceae           Tenacibaculum
## 8                 <NA>                    <NA>
## 131 Oceanospirillaceae          Neptuniibacter
## 306 Oceanospirillaceae             Marinomonas
## 186     Cryomorphaceae                Wandonia
## 29                <NA>                    <NA>
## 213  Flavobacteriaceae           Tenacibaculum
## 134 Oceanospirillaceae          Neptuniibacter
## 95      Cryomorphaceae             Owenweeksia
## 238  Flavobacteriaceae           Tenacibaculum
## 322   Alteromonadaceae         Aestuariibacter
## 393       Vibrionaceae                  Vibrio
## 391       Vibrionaceae                  Vibrio
## 343              SAR11 Candidatus Pelagibacter
## 227  Flavobacteriaceae           Tenacibaculum
## 297 Oceanospirillaceae             Marinomonas
## 15                <NA>                    <NA>
## 148 Oceanospirillaceae         Marinobacterium
## 237  Flavobacteriaceae           Tenacibaculum
## 370              SAR11 Candidatus Pelagibacter
## 217  Flavobacteriaceae           Tenacibaculum
## 394       Vibrionaceae                  Vibrio
## 209  Flavobacteriaceae           Tenacibaculum
## 313   Alteromonadaceae         Aestuariibacter
## 232  Flavobacteriaceae           Tenacibaculum
## 174     Cryomorphaceae                Wandonia
## 187     Cryomorphaceae                Wandonia
## 6                 <NA>                    <NA>
## 330   Alteromonadaceae         Aestuariibacter
## 392       Vibrionaceae                  Vibrio
## 303 Oceanospirillaceae             Marinomonas
## 386       Vibrionaceae                  Vibrio
## 319   Alteromonadaceae         Aestuariibacter
## 363              SAR11 Candidatus Pelagibacter
## 185     Cryomorphaceae                Wandonia
## 71      Cryomorphaceae             Owenweeksia
## 194     Cryomorphaceae                Wandonia
## 377       Vibrionaceae                  Vibrio
## 368              SAR11 Candidatus Pelagibacter
## 352              SAR11 Candidatus Pelagibacter
## 115 Oceanospirillaceae          Neptuniibacter
## 191     Cryomorphaceae                Wandonia
## 365              SAR11 Candidatus Pelagibacter
## 310   Alteromonadaceae         Aestuariibacter
## 221  Flavobacteriaceae           Tenacibaculum
## 373              SAR11 Candidatus Pelagibacter
## 280 Oceanospirillaceae             Marinomonas
## 334   Alteromonadaceae         Aestuariibacter
## 229  Flavobacteriaceae           Tenacibaculum
## 314   Alteromonadaceae         Aestuariibacter
## 146 Oceanospirillaceae         Marinobacterium
## 60   Flavobacteriaceae            Cellulophaga
## 324   Alteromonadaceae         Aestuariibacter
## 344              SAR11 Candidatus Pelagibacter
## 336   Alteromonadaceae         Aestuariibacter
## 11                <NA>                    <NA>
## 378       Vibrionaceae                  Vibrio
## 341              SAR11 Candidatus Pelagibacter
## 85      Cryomorphaceae             Owenweeksia
## 369              SAR11 Candidatus Pelagibacter
## 103 Oceanospirillaceae          Neptuniibacter
## 212  Flavobacteriaceae           Tenacibaculum
## 338   Alteromonadaceae         Aestuariibacter
## 248  Flavobacteriaceae              Aquibacter
## 308   Alteromonadaceae         Aestuariibacter
## 37   Flavobacteriaceae            Cellulophaga
## 192     Cryomorphaceae                Wandonia
## 169 Oceanospirillaceae         Marinobacterium
## 371              SAR11 Candidatus Pelagibacter
## 278 Oceanospirillaceae             Marinomonas
## 309   Alteromonadaceae         Aestuariibacter
## 298 Oceanospirillaceae             Marinomonas
## 302 Oceanospirillaceae             Marinomonas
## 318   Alteromonadaceae         Aestuariibacter
## 323   Alteromonadaceae         Aestuariibacter
## 132 Oceanospirillaceae          Neptuniibacter
## 23                <NA>                    <NA>
## 138 Oceanospirillaceae         Marinobacterium
## 197     Cryomorphaceae                Wandonia
## 135 Oceanospirillaceae          Neptuniibacter
## 307   Alteromonadaceae         Aestuariibacter
## 353              SAR11 Candidatus Pelagibacter
## 53   Flavobacteriaceae            Cellulophaga
## 77      Cryomorphaceae             Owenweeksia
## 355              SAR11 Candidatus Pelagibacter
## 364              SAR11 Candidatus Pelagibacter
## 335   Alteromonadaceae         Aestuariibacter
## 96      Cryomorphaceae             Owenweeksia
## 31                <NA>                    <NA>
## 38   Flavobacteriaceae            Cellulophaga
## 62   Flavobacteriaceae            Cellulophaga
## 36   Flavobacteriaceae            Cellulophaga
## 331   Alteromonadaceae         Aestuariibacter
## 46   Flavobacteriaceae            Cellulophaga
## 345              SAR11 Candidatus Pelagibacter
## 19                <NA>                    <NA>
## 315   Alteromonadaceae         Aestuariibacter
## 282 Oceanospirillaceae             Marinomonas
## 81      Cryomorphaceae             Owenweeksia
## 40   Flavobacteriaceae            Cellulophaga
## 74      Cryomorphaceae             Owenweeksia
## 113 Oceanospirillaceae          Neptuniibacter
## 275 Oceanospirillaceae             Marinomonas
## 167 Oceanospirillaceae         Marinobacterium
## 225  Flavobacteriaceae           Tenacibaculum
## 48   Flavobacteriaceae            Cellulophaga
## 243  Flavobacteriaceae              Aquibacter
## 159 Oceanospirillaceae         Marinobacterium
## 252  Flavobacteriaceae              Aquibacter
## 339   Alteromonadaceae         Aestuariibacter
## 78      Cryomorphaceae             Owenweeksia
## 173     Cryomorphaceae                Wandonia
## 88      Cryomorphaceae             Owenweeksia
## 140 Oceanospirillaceae         Marinobacterium
## 357              SAR11 Candidatus Pelagibacter
## 104 Oceanospirillaceae          Neptuniibacter
## 41   Flavobacteriaceae            Cellulophaga
## 154 Oceanospirillaceae         Marinobacterium
## 284 Oceanospirillaceae             Marinomonas
## 69      Cryomorphaceae             Owenweeksia
## 312   Alteromonadaceae         Aestuariibacter
## 202     Cryomorphaceae                Wandonia
## 172     Cryomorphaceae                Wandonia
## 72      Cryomorphaceae             Owenweeksia
## 359              SAR11 Candidatus Pelagibacter
## 164 Oceanospirillaceae         Marinobacterium
## 316   Alteromonadaceae         Aestuariibacter
## 325   Alteromonadaceae         Aestuariibacter
## 283 Oceanospirillaceae             Marinomonas
## 190     Cryomorphaceae                Wandonia
## 361              SAR11 Candidatus Pelagibacter
## 171     Cryomorphaceae                Wandonia
## 52   Flavobacteriaceae            Cellulophaga
## 198     Cryomorphaceae                Wandonia
## 139 Oceanospirillaceae         Marinobacterium
## 178     Cryomorphaceae                Wandonia
## 73      Cryomorphaceae             Owenweeksia
## 125 Oceanospirillaceae          Neptuniibacter
## 136 Oceanospirillaceae          Neptuniibacter
## 261  Flavobacteriaceae              Aquibacter
## 273 Oceanospirillaceae             Marinomonas
## 82      Cryomorphaceae             Owenweeksia
## 35   Flavobacteriaceae            Cellulophaga
## 123 Oceanospirillaceae          Neptuniibacter
## 79      Cryomorphaceae             Owenweeksia
## 76      Cryomorphaceae             Owenweeksia
## 245  Flavobacteriaceae              Aquibacter
## 372              SAR11 Candidatus Pelagibacter
## 83      Cryomorphaceae             Owenweeksia
## 264  Flavobacteriaceae              Aquibacter
## 137 Oceanospirillaceae         Marinobacterium
## 161 Oceanospirillaceae         Marinobacterium
## 204     Cryomorphaceae                Wandonia
## 347              SAR11 Candidatus Pelagibacter
## 63   Flavobacteriaceae            Cellulophaga
## 180     Cryomorphaceae                Wandonia
## 51   Flavobacteriaceae            Cellulophaga
## 265  Flavobacteriaceae              Aquibacter
## 242  Flavobacteriaceae              Aquibacter
## 260  Flavobacteriaceae              Aquibacter
## 327   Alteromonadaceae         Aestuariibacter
## 259  Flavobacteriaceae              Aquibacter
## 311   Alteromonadaceae         Aestuariibacter
## 75      Cryomorphaceae             Owenweeksia
## 241  Flavobacteriaceae              Aquibacter
## 367              SAR11 Candidatus Pelagibacter
## 346              SAR11 Candidatus Pelagibacter
## 244  Flavobacteriaceae              Aquibacter
## 294 Oceanospirillaceae             Marinomonas
## 249  Flavobacteriaceae              Aquibacter
## 152 Oceanospirillaceae         Marinobacterium
## 162 Oceanospirillaceae         Marinobacterium
## 240  Flavobacteriaceae              Aquibacter
## 145 Oceanospirillaceae         Marinobacterium
## 263  Flavobacteriaceae              Aquibacter
## 349              SAR11 Candidatus Pelagibacter
## 366              SAR11 Candidatus Pelagibacter
## 374              SAR11 Candidatus Pelagibacter
## 342              SAR11 Candidatus Pelagibacter
## 351              SAR11 Candidatus Pelagibacter
## 175     Cryomorphaceae                Wandonia
## 350              SAR11 Candidatus Pelagibacter
## 43   Flavobacteriaceae            Cellulophaga
## 200     Cryomorphaceae                Wandonia
## 250  Flavobacteriaceae              Aquibacter
## 360              SAR11 Candidatus Pelagibacter
## 133 Oceanospirillaceae          Neptuniibacter
## 124 Oceanospirillaceae          Neptuniibacter
## 356              SAR11 Candidatus Pelagibacter
## 47   Flavobacteriaceae            Cellulophaga
## 110 Oceanospirillaceae          Neptuniibacter
## 116 Oceanospirillaceae          Neptuniibacter
## 286 Oceanospirillaceae             Marinomonas
## 112 Oceanospirillaceae          Neptuniibacter
## 274 Oceanospirillaceae             Marinomonas
## 150 Oceanospirillaceae         Marinobacterium
## 196     Cryomorphaceae                Wandonia
## 55   Flavobacteriaceae            Cellulophaga
## 45   Flavobacteriaceae            Cellulophaga
## 156 Oceanospirillaceae         Marinobacterium
## 153 Oceanospirillaceae         Marinobacterium
## 160 Oceanospirillaceae         Marinobacterium
## 254  Flavobacteriaceae              Aquibacter
## 287 Oceanospirillaceae             Marinomonas
## 317   Alteromonadaceae         Aestuariibacter
## 87      Cryomorphaceae             Owenweeksia
## 234  Flavobacteriaceae           Tenacibaculum
## 42   Flavobacteriaceae            Cellulophaga
## 177     Cryomorphaceae                Wandonia
## 182     Cryomorphaceae                Wandonia
## 375       Vibrionaceae                  Vibrio
## 176     Cryomorphaceae                Wandonia
## 404       Vibrionaceae                  Vibrio
## 102     Cryomorphaceae             Owenweeksia
## 70      Cryomorphaceae             Owenweeksia
## 158 Oceanospirillaceae         Marinobacterium
## 383       Vibrionaceae                  Vibrio
## 328   Alteromonadaceae         Aestuariibacter
## 271  Flavobacteriaceae              Aquibacter
## 320   Alteromonadaceae         Aestuariibacter
## 385       Vibrionaceae                  Vibrio
## 239  Flavobacteriaceae              Aquibacter
## 108 Oceanospirillaceae          Neptuniibacter
## 258  Flavobacteriaceae              Aquibacter
## 381       Vibrionaceae                  Vibrio
## 157 Oceanospirillaceae         Marinobacterium
## 119 Oceanospirillaceae          Neptuniibacter
## 120 Oceanospirillaceae          Neptuniibacter
## 285 Oceanospirillaceae             Marinomonas
## 380       Vibrionaceae                  Vibrio
## 67   Flavobacteriaceae            Cellulophaga
## 80      Cryomorphaceae             Owenweeksia
## 109 Oceanospirillaceae          Neptuniibacter
## 118 Oceanospirillaceae          Neptuniibacter
## 170 Oceanospirillaceae         Marinobacterium
## 193     Cryomorphaceae                Wandonia
## 382       Vibrionaceae                  Vibrio
## 255  Flavobacteriaceae              Aquibacter
## 105 Oceanospirillaceae          Neptuniibacter
## 389       Vibrionaceae                  Vibrio
## 57   Flavobacteriaceae            Cellulophaga
## 304 Oceanospirillaceae             Marinomonas
## 58   Flavobacteriaceae            Cellulophaga
## 266  Flavobacteriaceae              Aquibacter
## 272  Flavobacteriaceae              Aquibacter
## 329   Alteromonadaceae         Aestuariibacter
## 49   Flavobacteriaceae            Cellulophaga
## 64   Flavobacteriaceae            Cellulophaga
## 65   Flavobacteriaceae            Cellulophaga
## 126 Oceanospirillaceae          Neptuniibacter
## 267  Flavobacteriaceae              Aquibacter
## 348              SAR11 Candidatus Pelagibacter
## 39   Flavobacteriaceae            Cellulophaga
## 44   Flavobacteriaceae            Cellulophaga
## 50   Flavobacteriaceae            Cellulophaga
## 54   Flavobacteriaceae            Cellulophaga
## 56   Flavobacteriaceae            Cellulophaga
## 59   Flavobacteriaceae            Cellulophaga
## 61   Flavobacteriaceae            Cellulophaga
## 66   Flavobacteriaceae            Cellulophaga
## 68   Flavobacteriaceae            Cellulophaga
## 89      Cryomorphaceae             Owenweeksia
## 90      Cryomorphaceae             Owenweeksia
## 91      Cryomorphaceae             Owenweeksia
## 92      Cryomorphaceae             Owenweeksia
## 93      Cryomorphaceae             Owenweeksia
## 94      Cryomorphaceae             Owenweeksia
## 97      Cryomorphaceae             Owenweeksia
## 100     Cryomorphaceae             Owenweeksia
## 101     Cryomorphaceae             Owenweeksia
## 106 Oceanospirillaceae          Neptuniibacter
## 107 Oceanospirillaceae          Neptuniibacter
## 114 Oceanospirillaceae          Neptuniibacter
## 121 Oceanospirillaceae          Neptuniibacter
## 122 Oceanospirillaceae          Neptuniibacter
## 127 Oceanospirillaceae          Neptuniibacter
## 128 Oceanospirillaceae          Neptuniibacter
## 129 Oceanospirillaceae          Neptuniibacter
## 130 Oceanospirillaceae          Neptuniibacter
## 142 Oceanospirillaceae         Marinobacterium
## 143 Oceanospirillaceae         Marinobacterium
## 144 Oceanospirillaceae         Marinobacterium
## 147 Oceanospirillaceae         Marinobacterium
## 149 Oceanospirillaceae         Marinobacterium
## 151 Oceanospirillaceae         Marinobacterium
## 155 Oceanospirillaceae         Marinobacterium
## 163 Oceanospirillaceae         Marinobacterium
## 165 Oceanospirillaceae         Marinobacterium
## 166 Oceanospirillaceae         Marinobacterium
## 179     Cryomorphaceae                Wandonia
## 183     Cryomorphaceae                Wandonia
## 184     Cryomorphaceae                Wandonia
## 188     Cryomorphaceae                Wandonia
## 189     Cryomorphaceae                Wandonia
## 195     Cryomorphaceae                Wandonia
## 201     Cryomorphaceae                Wandonia
## 246  Flavobacteriaceae              Aquibacter
## 247  Flavobacteriaceae              Aquibacter
## 251  Flavobacteriaceae              Aquibacter
## 253  Flavobacteriaceae              Aquibacter
## 257  Flavobacteriaceae              Aquibacter
## 262  Flavobacteriaceae              Aquibacter
## 268  Flavobacteriaceae              Aquibacter
## 269  Flavobacteriaceae              Aquibacter
## 270  Flavobacteriaceae              Aquibacter
## 291 Oceanospirillaceae             Marinomonas
## 300 Oceanospirillaceae             Marinomonas
## 301 Oceanospirillaceae             Marinomonas
## 321   Alteromonadaceae         Aestuariibacter
## 333   Alteromonadaceae         Aestuariibacter
## 340   Alteromonadaceae         Aestuariibacter
## 376       Vibrionaceae                  Vibrio
## 388       Vibrionaceae                  Vibrio
## 395       Vibrionaceae                  Vibrio
## 396       Vibrionaceae                  Vibrio
## 397       Vibrionaceae                  Vibrio
## 399       Vibrionaceae                  Vibrio
## 400       Vibrionaceae                  Vibrio
## 401       Vibrionaceae                  Vibrio
## 402       Vibrionaceae                  Vibrio
## 403       Vibrionaceae                  Vibrio
## 405       Vibrionaceae                  Vibrio
## 407       Vibrionaceae                  Vibrio
## 408       Vibrionaceae                  Vibrio
```

``` r
cg_filter1=subset_samples(cg_filt,innoc=="innoculum" | strain_26 =="A1" | strain_26==
                            "A2" | strain_26== "B1" | strain_26 =="C1")
physeq_merged <- merge_samples(cg_filter1, group = "definition_26")
```

```
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
```

``` r
glom=tax_glom(physeq_merged, "family")

abund=transform_sample_counts(glom, function(x) (x/sum(x)))
keep_taxa <- taxa_sums(abund) / nsamples(abund) > 0.01
ps_filtered <- prune_taxa(keep_taxa, abund)
ps_filtered=ps_prune(abund, min.samples = 0, min.reads = 0, min.abundance = 0.01)
```

```
## 82 features grouped as 'Others' in the output
```

``` r
tax_table(ps_filtered)
```

```
## Taxonomy Table:     [13 taxa by 7 taxonomic ranks]:
##        Kingdom    phyla                      class                        
## Otu6   "Bacteria" "Proteobacteria"           "Alphaproteobacteria"        
## Otu50  "Bacteria" "Proteobacteria"           "Alphaproteobacteria"        
## Otu68  "Bacteria" "Proteobacteria"           "Alphaproteobacteria"        
## Otu83  "Bacteria" "Bacteroidetes"            "Cytophagia"                 
## Otu2   "Bacteria" "Bacteroidetes"            "Flavobacteriia"             
## Otu26  "Bacteria" "candidate division WPS-2" "WPS-2_genera_incertae_sedis"
## Otu5   "Bacteria" "Proteobacteria"           "Gammaproteobacteria"        
## Otu7   "Bacteria" "Proteobacteria"           "Gammaproteobacteria"        
## Otu38  "Bacteria" "Proteobacteria"           "Betaproteobacteria"         
## Otu4   "Bacteria" "Proteobacteria"           "Gammaproteobacteria"        
## Otu36  "Bacteria" "Proteobacteria"           "Gammaproteobacteria"        
## Otu101 "Bacteria" "Proteobacteria"           "Gammaproteobacteria"        
## Others NA         NA                         NA                           
##        family                               genus species strain
## Otu6   "SAR11"                              NA    NA      NA    
## Otu50  "Rhodospirillales"                   NA    NA      NA    
## Otu68  "Rhodobacterales"                    NA    NA      NA    
## Otu83  "Cytophagales"                       NA    NA      NA    
## Otu2   "Flavobacteriales"                   NA    NA      NA    
## Otu26  "WPS-2_genera_incertae_sedis"        NA    NA      NA    
## Otu5   "Alteromonadales"                    NA    NA      NA    
## Otu7   "Vibrionales"                        NA    NA      NA    
## Otu38  "Methylophilales"                    NA    NA      NA    
## Otu4   "Oceanospirillales"                  NA    NA      NA    
## Otu36  "Pseudomonadales"                    NA    NA      NA    
## Otu101 "Gammaproteobacteria_incertae_sedis" NA    NA      NA    
## Others NA                                   NA    NA      NA
```

``` r
all_melt <- psmelt(ps_filtered)

all_melt <- all_melt %>%
  group_by(Sample, family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")
all_melt$family[is.na(all_melt$family)] <- "Other"
summary(as.factor(all_melt$family))
```

```
##                    Alteromonadales                       Cytophagales 
##                                  3                                  3 
##                   Flavobacteriales Gammaproteobacteria_incertae_sedis 
##                                  3                                  3 
##                    Methylophilales                  Oceanospirillales 
##                                  3                                  3 
##                              Other                    Pseudomonadales 
##                                  3                                  3 
##                    Rhodobacterales                   Rhodospirillales 
##                                  3                                  3 
##                              SAR11                        Vibrionales 
##                                  3                                  3 
##        WPS-2_genera_incertae_sedis 
##                                  3
```

``` r
length(summary(as.factor(all_melt$family)))
```

```
## [1] 13
```

``` r
all_melt <- all_melt %>%
  mutate(Sample = fct_relevel(Sample,
                              "Innoculum", "Coculture", "Algal-cell associated"))

all_figure=ggplot(all_melt, aes(x = Sample, y = Abundance, alluvium = family, stratum = family, fill = family)) +
  geom_flow(alpha = 0.7,na.rm=T) +
  geom_stratum() +
  theme_minimal(base_size=12) +
  theme(legend.position = 'none') +
  xlab("A1 sample flow") +
  ylab("Relative Abundance")
all_figure
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

``` r
library(ggalluvial)
library(ggplot2)


#A1 
cg_filter1=subset_samples(cg_filt,innoc=="innoculum" | strain_26 =="A1")
physeq_merged <- merge_samples(cg_filter1, group = "definition_26")
```

```
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
```

``` r
glom=tax_glom(physeq_merged, "family")
abund=transform_sample_counts(glom, function(x) (x/sum(x)))
a1_melt <- psmelt(abund)

keep_classes <- c("Alteromonadales", "Cytophagales", "Flavobacteriales",
                  "Gammaproteobacteria_incertae_sedis", "Methylophilales",
                  "Oceanospirillales", "Pseudomonadales", "Rhodobacterales",
                  "Rhodospirillales", "SAR11",
                  "Vibrionales", "WPS-2_genera_incertae_sedis")

a1_summed <- a1_melt %>%
  mutate(family = fct_other(family, keep = keep_classes, other_level = "Other"))

a1_summed <- a1_summed %>%
  group_by(Sample, family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")
library(forcats)
a1_summed <- a1_summed %>%
  mutate(Sample = fct_relevel(Sample,
                              "Innoculum", "Coculture", "Algal-cell associated"))
tol_muted <- c(
  "#88CCEE", "#44AA99", "#117733", "#332288",
  "#DDCC77", "#999933", "#CC6677", "#882255",
  "#AA4499", "#DDDDDD", "#6699CC", "#888888", 'black'
)

a1_figure=ggplot(a1_summed, aes(x = Sample, y = Abundance, alluvium = family, stratum = family, fill = family)) +
  geom_flow(alpha = 0.7,na.rm=T) +
  geom_stratum() +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'none') +
  xlab("A1") +
  scale_fill_manual(values = tol_muted) +
  ylab("Relative Abundance")
a1_figure
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-16-2.png)<!-- -->

``` r
a1_figure2=ggplot(a1_summed, aes(x = Sample, y = Abundance, alluvium = family, stratum = family, fill = family)) +
  geom_flow(alpha = 0.7,na.rm=T) +
  geom_stratum() +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  xlab("A1") +
  scale_fill_manual(values = tol_muted) +
  ylab("Relative Abundance")
a1_figure2
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-16-3.png)<!-- -->

``` r
#ggsave('flowplot_legend.svg', plot=a1_figure2)

#A2
cg_filter1=subset_samples(cg_filt,innoc=="innoculum" | strain_26 =="A2")
physeq_merged <- merge_samples(cg_filter1, group = "definition_26")
```

```
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
```

``` r
glom=tax_glom(physeq_merged, "family")
abund=transform_sample_counts(glom, function(x) (x/sum(x)))
a2_melt <- psmelt(abund)

a2_summed <- a2_melt %>%
  mutate(family = fct_other(family, keep = keep_classes, other_level = "Other"))

a2_summed <- a2_summed %>%
  group_by(Sample, family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

library(forcats)
a2_summed <- a2_summed %>%
  mutate(Sample = fct_relevel(Sample,
                                     "Innoculum", "Coculture", "Algal-cell associated"))

a2_figure=ggplot(a2_summed, aes(x = Sample, y = Abundance, alluvium = family, stratum = family, fill = family)) +
  geom_flow(alpha = 0.7) +
  geom_stratum() +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'none') +
  xlab("A2") +
  scale_fill_manual(values = tol_muted) +
  ylab("Relative Abundance")
a2_figure
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-16-4.png)<!-- -->

``` r
#B1
cg_filter1=subset_samples(cg_filt,innoc=="innoculum" | strain_26 =="B1")
physeq_merged <- merge_samples(cg_filter1, group = "definition_26")
```

```
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
```

``` r
glom=tax_glom(physeq_merged, "family")
abund=transform_sample_counts(glom, function(x) (x/sum(x)))
b1_melt <- psmelt(abund)
b1_summed <- b1_melt %>%
  mutate(family = fct_other(family, keep = keep_classes, other_level = "Other"))

b1_summed <- b1_summed %>%
  group_by(Sample, family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

b1_summed <- b1_summed %>%
  mutate(Sample = fct_relevel(Sample,
                              "Innoculum", "Coculture", "Algal-cell associated"))

b1_figure=ggplot(b1_summed, aes(x = Sample, y = Abundance, alluvium = family, stratum = family, fill = family)) +
  geom_flow(alpha = 0.7) +
  geom_stratum() +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'none') +
  xlab("B1") +
  scale_fill_manual(values = tol_muted) +
  ylab("Relative Abundance")
b1_figure
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-16-5.png)<!-- -->

``` r
# C1
cg_filter1=subset_samples(cg_filt,innoc=="innoculum" | strain_26 =="C1")
physeq_merged <- merge_samples(cg_filter1, group = "definition_26")
```

```
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
```

``` r
glom=tax_glom(physeq_merged, "family")
abund=transform_sample_counts(glom, function(x) (x/sum(x)))

c1_melt <- psmelt(abund)
c1_summed <- c1_melt %>%
  mutate(family = fct_other(family, keep = keep_classes, other_level = "Other"))
c1_summed <- c1_summed %>%
  group_by(Sample, family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

c1_summed <- c1_summed %>%
  mutate(Sample = fct_relevel(Sample,
                              "Innoculum", "Coculture", "Algal-cell associated"))

c1_figure=ggplot(c1_summed, aes(x = Sample, y = Abundance, alluvium = family, stratum = family, fill = family)) +
  geom_flow(alpha = 0.7) +
  geom_stratum() +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'none') +
  xlab("C1") +
  scale_fill_manual(values = tol_muted) +
  ylab("Relative Abundance")
c1_figure
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-16-6.png)<!-- -->

``` r
cowplot::plot_grid(a1_figure, a2_figure, b1_figure,c1_figure)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-16-7.png)<!-- -->

``` r
cowww=cowplot::plot_grid(a1_figure, a2_figure, b1_figure,c1_figure,ncol =4)

#ggsave('ribbon_plot.svg', plot=cowww, dpi=300, height=4, width=16)
#cowplot::plot_grid(a1_figure, a2_figure, b1_figure,c1_figure, 
 #                  ncol=1)
```




# ribbon plot genus


``` r
library(ggalluvial)
library(ggplot2)
df_melt
```

```
##        OTU Sample    Abundance X.OTU.ID      info strain_26 diatom_strain
## 296   Otu4  WGA18 0.6501439737    WGA18 WGA_RD4-A        A2        NbO-D4
## 295   Otu4  WGA19 0.6147675854    WGA19 WGA_RD4-A        A2        NbO-D4
## 281   Otu4  OMA22 0.5571781160    OMA22  RD4-C_18        A2        NbO-D4
## 111  Otu12  WGA33 0.5277663513    WGA33 WGA_NB7-A        C1        NbQ-B7
## 18  Others  WGA40 0.5168654875    WGA40 WGA_NB7-B        C1        NbQ-B7
## 224   Otu2  WGA34 0.5143973673    WGA34 WGA_NB7-A        C1        NbQ-B7
## 231   Otu2  OMA20 0.5139860140    OMA20  RD4-A_18        A2        NbO-D4
## 222   Otu2  OMA13 0.5111065405    OMA13  NB7-C_18        C1        NbQ-B7
## 7   Others  WGA39 0.5074043603    WGA39 WGA_NB7-B        C1        NbQ-B7
## 279   Otu4  WGA20 0.4845742493    WGA20 WGA_RD4-A        A2        NbO-D4
## 214   Otu2  OMA12 0.4790209790    OMA12  NB7-B_18        C1        NbQ-B7
## 358   Otu6  WGA41 0.4625668449    WGA41 WGA_NB7-B        C1        NbQ-B7
## 390   Otu7  OMA19 0.4574249280    OMA19  RA5-C_18        A1        NbO-A5
## 205   Otu2  OMA26 0.4477581242    OMA26  YE5-A_18        B1       NbP-YE5
## 207   Otu2  OMA27 0.4224598930    OMA27  YE5-B_18        B1       NbP-YE5
## 289   Otu4  WGA53 0.3879062114    WGA53 WGA_YE5-B        B1       NbP-YE5
## 362   Otu6  WGA54 0.3835870012    WGA54 WGA_YE5-B        B1       NbP-YE5
## 206   Otu2  OMA11 0.3745372275    OMA11  NB7-A_18        C1        NbQ-B7
## 406   Otu7   WGA4 0.3712464007     WGA4 WGA_RA5-A        A1        NbO-A5
## 299   Otu4  WGA10 0.3679555738    WGA10 WGA_RA5-B        A1        NbO-A5
## 12  Others  OMA28 0.3455368161    OMA28  YE5-C_18        B1       NbP-YE5
## 290   Otu4  WGA27 0.3443027561    WGA27 WGA_RD4-B        A2        NbO-D4
## 379   Otu7  OMA17 0.3414232826    OMA17  RA5-A_18        A1        NbO-A5
## 211   Otu2  OMA28 0.3364870424    OMA28  YE5-C_18        B1       NbP-YE5
## 34  Others  WGA55 0.3362813657    WGA55 WGA_YE5-B        B1       NbP-YE5
## 141  Otu14  WGA25 0.3354586590    WGA25 WGA_RD4-B        A2        NbO-D4
## 220   Otu2  OMA18 0.3280542986    OMA18  RA5-B_18        A1        NbO-A5
## 210   Otu2  OMA21 0.3276429453    OMA21  RD4-B_18        A2        NbO-D4
## 99   Otu11  WGA34 0.3175647882    WGA34 WGA_NB7-A        C1        NbQ-B7
## 326   Otu5  WGA46 0.3171534348    WGA46 WGA_YE5-A        B1       NbP-YE5
## 384   Otu7   WGA3 0.3140682847     WGA3 WGA_RA5-A        A1        NbO-A5
## 292   Otu4  WGA46 0.3138626080    WGA46 WGA_YE5-A        B1       NbP-YE5
## 5   Others  OMA11 0.3027560675    OMA11  NB7-A_18        C1        NbQ-B7
## 13  Others  OMA12 0.2982311806    OMA12  NB7-B_18        C1        NbQ-B7
## 256  Otu24  WGA48 0.2914438503    WGA48 WGA_YE5-A        B1       NbP-YE5
## 305   Otu4  WGA55 0.2815713698    WGA55 WGA_YE5-B        B1       NbP-YE5
## 181  Otu19  WGA45 0.2782805430    WGA45 WGA_YE5-A        B1       NbP-YE5
## 277   Otu4  OMA21 0.2747840395    OMA21  RD4-B_18        A2        NbO-D4
## 28  Others  OMA26 0.2556561086    OMA26  YE5-A_18        B1       NbP-YE5
## 1   Others   WGA5 0.2542163719     WGA5 WGA_RA5-A        A1        NbO-A5
## 2   Others  WGA45 0.2519539284    WGA45 WGA_YE5-A        B1       NbP-YE5
## 30  Others  WGA41 0.2441382147    WGA41 WGA_NB7-B        C1        NbQ-B7
## 32  Others  OMA20 0.2414644179    OMA20  RD4-A_18        A2        NbO-D4
## 168  Otu14  WGA41 0.2400246812    WGA41 WGA_NB7-B        C1        NbQ-B7
## 21  Others  OMA18 0.2398190045    OMA18  RA5-B_18        A1        NbO-A5
## 235   Otu2  WGA12 0.2381735911    WGA12 WGA_RA5-B        A1        NbO-A5
## 4   Others  OMA27 0.2354997943    OMA27  YE5-B_18        B1       NbP-YE5
## 10  Others  WGA53 0.2334430276    WGA53 WGA_YE5-B        B1       NbP-YE5
## 86   Otu11  WGA12 0.2305635541    WGA12 WGA_RA5-B        A1        NbO-A5
## 25  Others  OMA13 0.2295351707    OMA13  NB7-C_18        C1        NbQ-B7
## 20  Others  WGA13 0.2295351707    WGA13 WGA_RA5-B        A1        NbO-A5
## 387   Otu7  WGA25 0.2266556972    WGA25 WGA_RD4-B        A2        NbO-D4
## 208   Otu2  OMA17 0.2225421637    OMA17  RA5-A_18        A1        NbO-A5
## 218   Otu2  WGA27 0.2219251337    WGA27 WGA_RD4-B        A2        NbO-D4
## 236   Otu2  WGA39 0.2219251337    WGA39 WGA_NB7-B        C1        NbQ-B7
## 14  Others  WGA46 0.2186343069    WGA46 WGA_YE5-A        B1       NbP-YE5
## 26  Others  WGA48 0.2159605101    WGA48 WGA_YE5-A        B1       NbP-YE5
## 219   Otu2  OMA19 0.2108185932    OMA19  RA5-C_18        A1        NbO-A5
## 398   Otu7   WGA5 0.2104072398     WGA5 WGA_RA5-A        A1        NbO-A5
## 9   Others  OMA17 0.2001234060    OMA17  RA5-A_18        A1        NbO-A5
## 288   Otu4  WGA33 0.1995063760    WGA33 WGA_NB7-A        C1        NbQ-B7
## 24  Others  WGA10 0.1879884821    WGA10 WGA_RA5-B        A1        NbO-A5
## 276   Otu4   WGA3 0.1863430687     WGA3 WGA_RA5-A        A1        NbO-A5
## 228   Otu2  WGA48 0.1803784451    WGA48 WGA_YE5-A        B1       NbP-YE5
## 223   Otu2  WGA10 0.1793500617    WGA10 WGA_RA5-B        A1        NbO-A5
## 22  Others  WGA54 0.1762649116    WGA54 WGA_YE5-B        B1       NbP-YE5
## 3   Others   WGA3 0.1737967914     WGA3 WGA_RA5-A        A1        NbO-A5
## 199  Otu19  WGA19 0.1723570547    WGA19 WGA_RD4-A        A2        NbO-D4
## 332   Otu5  WGA12 0.1713286713    WGA12 WGA_RA5-B        A1        NbO-A5
## 17  Others  OMA21 0.1709173180    OMA21  RD4-B_18        A2        NbO-D4
## 16  Others  WGA20 0.1707116413    WGA20 WGA_RD4-A        A2        NbO-D4
## 216   Otu2  WGA20 0.1655697244    WGA20 WGA_RD4-A        A2        NbO-D4
## 337   Otu5   WGA5 0.1631016043     WGA5 WGA_RA5-A        A1        NbO-A5
## 230   Otu2  WGA40 0.1606334842    WGA40 WGA_NB7-B        C1        NbQ-B7
## 33  Others  OMA19 0.1585767174    OMA19  RA5-C_18        A1        NbO-A5
## 226   Otu2  WGA53 0.1559029206    WGA53 WGA_YE5-B        B1       NbP-YE5
## 203  Otu19  WGA13 0.1499382970    WGA13 WGA_RA5-B        A1        NbO-A5
## 354   Otu6  WGA48 0.1497326203    WGA48 WGA_YE5-A        B1       NbP-YE5
## 84   Otu11  WGA18 0.1495269436    WGA18 WGA_RD4-A        A2        NbO-D4
## 117  Otu12  WGA45 0.1472645002    WGA45 WGA_YE5-A        B1       NbP-YE5
## 215   Otu2  WGA25 0.1470588235    WGA25 WGA_RD4-B        A2        NbO-D4
## 27  Others  WGA27 0.1460304401    WGA27 WGA_RD4-B        A2        NbO-D4
## 293   Otu4  WGA13 0.1431509667    WGA13 WGA_RA5-B        A1        NbO-A5
## 98   Otu11  WGA40 0.1404771699    WGA40 WGA_NB7-B        C1        NbQ-B7
## 233   Otu2  WGA18 0.1402714932    WGA18 WGA_RD4-A        A2        NbO-D4
## 8   Others  WGA12 0.1394487865    WGA12 WGA_RA5-B        A1        NbO-A5
## 131  Otu12  WGA39 0.1375976964    WGA39 WGA_NB7-B        C1        NbQ-B7
## 306   Otu4   WGA4 0.1367749897     WGA4 WGA_RA5-A        A1        NbO-A5
## 186  Otu19  OMA18 0.1330728095    OMA18  RA5-B_18        A1        NbO-A5
## 29  Others  OMA22 0.1330728095    OMA22  RD4-C_18        A2        NbO-D4
## 213   Otu2   WGA3 0.1314273961     WGA3 WGA_RA5-A        A1        NbO-A5
## 134  Otu12  WGA53 0.1312217195    WGA53 WGA_YE5-B        B1       NbP-YE5
## 95   Otu11   WGA4 0.1287535993     WGA4 WGA_RA5-A        A1        NbO-A5
## 238   Otu2   WGA4 0.1256684492     WGA4 WGA_RA5-A        A1        NbO-A5
## 322   Otu5  WGA54 0.1248457425    WGA54 WGA_YE5-B        B1       NbP-YE5
## 393   Otu7  WGA12 0.1246400658    WGA12 WGA_RA5-B        A1        NbO-A5
## 391   Otu7  OMA18 0.1229946524    OMA18  RA5-B_18        A1        NbO-A5
## 343   Otu6  WGA45 0.1225832991    WGA45 WGA_YE5-A        B1       NbP-YE5
## 227   Otu2   WGA5 0.1223776224     WGA5 WGA_RA5-A        A1        NbO-A5
## 297   Otu4  WGA54 0.1209378856    WGA54 WGA_YE5-B        B1       NbP-YE5
## 15  Others  WGA25 0.1194981489    WGA25 WGA_RD4-B        A2        NbO-D4
## 148  Otu14  WGA13 0.1180584122    WGA13 WGA_RA5-B        A1        NbO-A5
## 237   Otu2  WGA55 0.1164129988    WGA55 WGA_YE5-B        B1       NbP-YE5
## 370   Otu6  WGA55 0.1129164953    WGA55 WGA_YE5-B        B1       NbP-YE5
## 217   Otu2  WGA45 0.1125051419    WGA45 WGA_YE5-A        B1       NbP-YE5
## 394   Otu7  WGA10 0.1120937886    WGA10 WGA_RA5-B        A1        NbO-A5
## 209   Otu2  OMA22 0.1116824352    OMA22  RD4-C_18        A2        NbO-D4
## 313   Otu5  WGA33 0.1104483752    WGA33 WGA_NB7-A        C1        NbQ-B7
## 232   Otu2  WGA19 0.1071575483    WGA19 WGA_RD4-A        A2        NbO-D4
## 174  Otu19  OMA17 0.1057178116    OMA17  RA5-A_18        A1        NbO-A5
## 187  Otu19  WGA54 0.1030440148    WGA54 WGA_YE5-B        B1       NbP-YE5
## 6   Others   WGA4 0.1013986014     WGA4 WGA_RA5-A        A1        NbO-A5
## 330   Otu5  WGA13 0.1001645413    WGA13 WGA_RA5-B        A1        NbO-A5
## 392   Otu7  WGA13 0.0993418346    WGA13 WGA_RA5-B        A1        NbO-A5
## 303   Otu4   WGA5 0.0964623612     WGA5 WGA_RA5-A        A1        NbO-A5
## 386   Otu7  WGA27 0.0917317976    WGA27 WGA_RD4-B        A2        NbO-D4
## 319   Otu5  OMA13 0.0878239408    OMA13  NB7-C_18        C1        NbQ-B7
## 363   Otu6  WGA34 0.0876182641    WGA34 WGA_NB7-A        C1        NbQ-B7
## 185  Otu19   WGA4 0.0841217606     WGA4 WGA_RA5-A        A1        NbO-A5
## 71   Otu11  OMA18 0.0839160839    OMA18  RA5-B_18        A1        NbO-A5
## 194  Otu19  WGA25 0.0830933772    WGA25 WGA_RD4-B        A2        NbO-D4
## 377   Otu7  OMA22 0.0818593172    OMA22  RD4-C_18        A2        NbO-D4
## 368   Otu6  OMA27 0.0812422871    OMA27  YE5-B_18        B1       NbP-YE5
## 352   Otu6  OMA28 0.0793911970    OMA28  YE5-C_18        B1       NbP-YE5
## 115  Otu12  WGA48 0.0781571370    WGA48 WGA_YE5-A        B1       NbP-YE5
## 191  Otu19   WGA5 0.0765117236     WGA5 WGA_RA5-A        A1        NbO-A5
## 365   Otu6  OMA13 0.0765117236    OMA13  NB7-C_18        C1        NbQ-B7
## 310   Otu5  OMA28 0.0758946935    OMA28  YE5-C_18        B1       NbP-YE5
## 221   Otu2  WGA13 0.0752776635    WGA13 WGA_RA5-B        A1        NbO-A5
## 373   Otu6  OMA20 0.0744549568    OMA20  RD4-A_18        A2        NbO-D4
## 280   Otu4  WGA45 0.0738379268    WGA45 WGA_YE5-A        B1       NbP-YE5
## 334   Otu5  WGA55 0.0736322501    WGA55 WGA_YE5-B        B1       NbP-YE5
## 229   Otu2  WGA46 0.0734265734    WGA46 WGA_YE5-A        B1       NbP-YE5
## 314   Otu5  OMA26 0.0709584533    OMA26  YE5-A_18        B1       NbP-YE5
## 146  Otu14  OMA11 0.0707527766    OMA11  NB7-A_18        C1        NbQ-B7
## 60   Otu10  WGA39 0.0699300699    WGA39 WGA_NB7-B        C1        NbQ-B7
## 324   Otu5  WGA27 0.0695187166    WGA27 WGA_RD4-B        A2        NbO-D4
## 344   Otu6  WGA20 0.0691073632    WGA20 WGA_RD4-A        A2        NbO-D4
## 336   Otu5  WGA53 0.0689016865    WGA53 WGA_YE5-B        B1       NbP-YE5
## 11  Others  WGA33 0.0682846565    WGA33 WGA_NB7-A        C1        NbQ-B7
## 378   Otu7  OMA21 0.0674619498    OMA21  RD4-B_18        A2        NbO-D4
## 341   Otu6  OMA11 0.0658165364    OMA11  NB7-A_18        C1        NbQ-B7
## 85   Otu11  WGA13 0.0656108597    WGA13 WGA_RA5-B        A1        NbO-A5
## 369   Otu6  OMA26 0.0656108597    OMA26  YE5-A_18        B1       NbP-YE5
## 103  Otu12  OMA27 0.0651995064    OMA27  YE5-B_18        B1       NbP-YE5
## 212   Otu2  WGA33 0.0651995064    WGA33 WGA_NB7-A        C1        NbQ-B7
## 338   Otu5  WGA48 0.0631427396    WGA48 WGA_YE5-A        B1       NbP-YE5
## 248  Otu24  WGA20 0.0625257096    WGA20 WGA_RD4-A        A2        NbO-D4
## 308   Otu5  OMA27 0.0625257096    OMA27  YE5-B_18        B1       NbP-YE5
## 37   Otu10  OMA26 0.0617030029    OMA26  YE5-A_18        B1       NbP-YE5
## 192  Otu19   WGA3 0.0573837927     WGA3 WGA_RA5-A        A1        NbO-A5
## 169  Otu14  WGA40 0.0557383793    WGA40 WGA_NB7-B        C1        NbQ-B7
## 371   Otu6  WGA39 0.0555327026    WGA39 WGA_NB7-B        C1        NbQ-B7
## 278   Otu4  WGA25 0.0547099959    WGA25 WGA_RD4-B        A2        NbO-D4
## 309   Otu5  OMA12 0.0547099959    OMA12  NB7-B_18        C1        NbQ-B7
## 298   Otu4  WGA12 0.0545043192    WGA12 WGA_RA5-B        A1        NbO-A5
## 302   Otu4  WGA40 0.0530645825    WGA40 WGA_NB7-B        C1        NbQ-B7
## 318   Otu5  OMA20 0.0526532291    OMA20  RD4-A_18        A2        NbO-D4
## 323   Otu5   WGA3 0.0514191691     WGA3 WGA_RA5-A        A1        NbO-A5
## 132  Otu12  WGA55 0.0501851090    WGA55 WGA_YE5-B        B1       NbP-YE5
## 23  Others  WGA34 0.0495680790    WGA34 WGA_NB7-A        C1        NbQ-B7
## 138  Otu14  WGA27 0.0495680790    WGA27 WGA_RD4-B        A2        NbO-D4
## 197  Otu19  OMA20 0.0485396956    OMA20  RD4-A_18        A2        NbO-D4
## 135  Otu12  WGA41 0.0481283422    WGA41 WGA_NB7-B        C1        NbQ-B7
## 307   Otu5  OMA11 0.0481283422    OMA11  NB7-A_18        C1        NbQ-B7
## 353   Otu6  OMA12 0.0479226656    OMA12  NB7-B_18        C1        NbQ-B7
## 53   Otu10  WGA46 0.0468942822    WGA46 WGA_YE5-A        B1       NbP-YE5
## 77   Otu11   WGA3 0.0462772522     WGA3 WGA_RA5-A        A1        NbO-A5
## 355   Otu6  WGA27 0.0450431921    WGA27 WGA_RD4-B        A2        NbO-D4
## 364   Otu6  WGA10 0.0444261621    WGA10 WGA_RA5-B        A1        NbO-A5
## 335   Otu5  WGA40 0.0431921020    WGA40 WGA_NB7-B        C1        NbQ-B7
## 96   Otu11  WGA10 0.0429864253    WGA10 WGA_RA5-B        A1        NbO-A5
## 31  Others  WGA19 0.0419580420    WGA19 WGA_RD4-A        A2        NbO-D4
## 38   Otu10  OMA12 0.0409296586    OMA12  NB7-B_18        C1        NbQ-B7
## 62   Otu10  WGA54 0.0407239819    WGA54 WGA_YE5-B        B1       NbP-YE5
## 36   Otu10  OMA11 0.0403126285    OMA11  NB7-A_18        C1        NbQ-B7
## 331   Otu5  OMA17 0.0399012752    OMA17  RA5-A_18        A1        NbO-A5
## 46   Otu10  OMA13 0.0394899218    OMA13  NB7-C_18        C1        NbQ-B7
## 345   Otu6  OMA21 0.0394899218    OMA21  RD4-B_18        A2        NbO-D4
## 19  Others  WGA18 0.0392842452    WGA18 WGA_RD4-A        A2        NbO-D4
## 315   Otu5  OMA19 0.0392842452    OMA19  RA5-C_18        A1        NbO-A5
## 282   Otu4  OMA19 0.0390785685    OMA19  RA5-C_18        A1        NbO-A5
## 81   Otu11  WGA20 0.0372274784    WGA20 WGA_RD4-A        A2        NbO-D4
## 40   Otu10  OMA27 0.0364047717    OMA27  YE5-B_18        B1       NbP-YE5
## 74   Otu11  OMA22 0.0361990950    OMA22  RD4-C_18        A2        NbO-D4
## 113  Otu12  OMA11 0.0359934183    OMA11  NB7-A_18        C1        NbQ-B7
## 275   Otu4  OMA27 0.0355820650    OMA27  YE5-B_18        B1       NbP-YE5
## 167  Otu14   WGA5 0.0347593583     WGA5 WGA_RA5-A        A1        NbO-A5
## 225   Otu2  WGA54 0.0347593583    WGA54 WGA_YE5-B        B1       NbP-YE5
## 48   Otu10  OMA28 0.0343480049    OMA28  YE5-C_18        B1       NbP-YE5
## 243  Otu24  OMA11 0.0339366516    OMA11  NB7-A_18        C1        NbQ-B7
## 159  Otu14  WGA19 0.0327025915    WGA19 WGA_RD4-A        A2        NbO-D4
## 252  Otu24  OMA28 0.0320855615    OMA28  YE5-C_18        B1       NbP-YE5
## 339   Otu5   WGA4 0.0320855615     WGA4 WGA_RA5-A        A1        NbO-A5
## 78   Otu11   WGA5 0.0312628548     WGA5 WGA_RA5-A        A1        NbO-A5
## 173  Otu19  OMA21 0.0312628548    OMA21  RD4-B_18        A2        NbO-D4
## 88   Otu11  OMA28 0.0296174414    OMA28  YE5-C_18        B1       NbP-YE5
## 140  Otu14  OMA28 0.0294117647    OMA28  YE5-C_18        B1       NbP-YE5
## 357   Otu6  OMA22 0.0294117647    OMA22  RD4-C_18        A2        NbO-D4
## 104  Otu12  OMA26 0.0290004114    OMA26  YE5-A_18        B1       NbP-YE5
## 41   Otu10  OMA20 0.0287947347    OMA20  RD4-A_18        A2        NbO-D4
## 154  Otu14  OMA12 0.0285890580    OMA12  NB7-B_18        C1        NbQ-B7
## 284   Otu4  OMA17 0.0285890580    OMA17  RA5-A_18        A1        NbO-A5
## 69   Otu11  OMA21 0.0281777046    OMA21  RD4-B_18        A2        NbO-D4
## 312   Otu5  OMA21 0.0279720280    OMA21  RD4-B_18        A2        NbO-D4
## 202  Otu19  WGA10 0.0277663513    WGA10 WGA_RA5-B        A1        NbO-A5
## 172  Otu19  OMA22 0.0275606746    OMA22  RD4-C_18        A2        NbO-D4
## 72   Otu11  OMA17 0.0271493213    OMA17  RA5-A_18        A1        NbO-A5
## 359   Otu6  WGA19 0.0263266146    WGA19 WGA_RD4-A        A2        NbO-D4
## 164  Otu14  WGA55 0.0248868778    WGA55 WGA_YE5-B        B1       NbP-YE5
## 316   Otu5  OMA18 0.0246812012    OMA18  RA5-B_18        A1        NbO-A5
## 325   Otu5  WGA25 0.0246812012    WGA25 WGA_RD4-B        A2        NbO-D4
## 283   Otu4  OMA18 0.0242698478    OMA18  RA5-B_18        A1        NbO-A5
## 190  Otu19  OMA28 0.0230357877    OMA28  YE5-C_18        B1       NbP-YE5
## 361   Otu6  OMA19 0.0230357877    OMA19  RA5-C_18        A1        NbO-A5
## 171  Otu19  OMA26 0.0220074044    OMA26  YE5-A_18        B1       NbP-YE5
## 52   Otu10  WGA48 0.0211846977    WGA48 WGA_YE5-A        B1       NbP-YE5
## 198  Otu19  OMA19 0.0203619910    OMA19  RA5-C_18        A1        NbO-A5
## 139  Otu14  WGA33 0.0199506376    WGA33 WGA_NB7-A        C1        NbQ-B7
## 178  Otu19  OMA27 0.0199506376    OMA27  YE5-B_18        B1       NbP-YE5
## 73   Otu11  OMA13 0.0191279309    OMA13  NB7-C_18        C1        NbQ-B7
## 125  Otu12  OMA12 0.0180995475    OMA12  NB7-B_18        C1        NbQ-B7
## 136  Otu12  WGA40 0.0180995475    WGA40 WGA_NB7-B        C1        NbQ-B7
## 261  Otu24  WGA13 0.0180995475    WGA13 WGA_RA5-B        A1        NbO-A5
## 273   Otu4  OMA26 0.0176881942    OMA26  YE5-A_18        B1       NbP-YE5
## 82   Otu11  OMA20 0.0174825175    OMA20  RD4-A_18        A2        NbO-D4
## 35   Otu10  OMA21 0.0172768408    OMA21  RD4-B_18        A2        NbO-D4
## 123  Otu12  WGA34 0.0170711641    WGA34 WGA_NB7-A        C1        NbQ-B7
## 79   Otu11  WGA27 0.0168654875    WGA27 WGA_RD4-B        A2        NbO-D4
## 76   Otu11  OMA27 0.0166598108    OMA27  YE5-B_18        B1       NbP-YE5
## 245  Otu24  OMA26 0.0166598108    OMA26  YE5-A_18        B1       NbP-YE5
## 372   Otu6  WGA12 0.0166598108    WGA12 WGA_RA5-B        A1        NbO-A5
## 83   Otu11  OMA19 0.0160427807    OMA19  RA5-C_18        A1        NbO-A5
## 264  Otu24  WGA10 0.0160427807    WGA10 WGA_RA5-B        A1        NbO-A5
## 137  Otu14   WGA3 0.0152200740     WGA3 WGA_RA5-A        A1        NbO-A5
## 161  Otu14  OMA19 0.0152200740    OMA19  RA5-C_18        A1        NbO-A5
## 204  Otu19  WGA12 0.0152200740    WGA12 WGA_RA5-B        A1        NbO-A5
## 347   Otu6   WGA4 0.0150143974     WGA4 WGA_RA5-A        A1        NbO-A5
## 63   Otu10  WGA53 0.0143973673    WGA53 WGA_YE5-B        B1       NbP-YE5
## 180  Otu19  WGA46 0.0139860140    WGA46 WGA_YE5-A        B1       NbP-YE5
## 51   Otu10  WGA27 0.0137803373    WGA27 WGA_RD4-B        A2        NbO-D4
## 265  Otu24  WGA34 0.0137803373    WGA34 WGA_NB7-A        C1        NbQ-B7
## 242  Otu24  OMA12 0.0133689840    OMA12  NB7-B_18        C1        NbQ-B7
## 260  Otu24  OMA18 0.0133689840    OMA18  RA5-B_18        A1        NbO-A5
## 327   Otu5  WGA45 0.0133689840    WGA45 WGA_YE5-A        B1       NbP-YE5
## 259  Otu24  OMA19 0.0125462773    OMA19  RA5-C_18        A1        NbO-A5
## 311   Otu5  OMA22 0.0119292472    OMA22  RD4-C_18        A2        NbO-D4
## 75   Otu11  OMA11 0.0117235705    OMA11  NB7-A_18        C1        NbQ-B7
## 241  Otu24  OMA13 0.0117235705    OMA13  NB7-C_18        C1        NbQ-B7
## 367   Otu6   WGA3 0.0117235705     WGA3 WGA_RA5-A        A1        NbO-A5
## 346   Otu6  WGA40 0.0115178939    WGA40 WGA_NB7-B        C1        NbQ-B7
## 244  Otu24  OMA27 0.0113122172    OMA27  YE5-B_18        B1       NbP-YE5
## 294   Otu4  OMA20 0.0113122172    OMA20  RD4-A_18        A2        NbO-D4
## 249  Otu24  OMA17 0.0111065405    OMA17  RA5-A_18        A1        NbO-A5
## 152  Otu14  WGA10 0.0106951872    WGA10 WGA_RA5-B        A1        NbO-A5
## 162  Otu14  OMA18 0.0104895105    OMA18  RA5-B_18        A1        NbO-A5
## 240  Otu24  OMA21 0.0104895105    OMA21  RD4-B_18        A2        NbO-D4
## 145  Otu14  OMA27 0.0094611271    OMA27  YE5-B_18        B1       NbP-YE5
## 263  Otu24  WGA12 0.0094611271    WGA12 WGA_RA5-B        A1        NbO-A5
## 349   Otu6  OMA18 0.0094611271    OMA18  RA5-B_18        A1        NbO-A5
## 366   Otu6   WGA5 0.0094611271     WGA5 WGA_RA5-A        A1        NbO-A5
## 374   Otu6  OMA17 0.0092554504    OMA17  RA5-A_18        A1        NbO-A5
## 342   Otu6  WGA46 0.0086384204    WGA46 WGA_YE5-A        B1       NbP-YE5
## 351   Otu6  WGA33 0.0086384204    WGA33 WGA_NB7-A        C1        NbQ-B7
## 175  Otu19  OMA13 0.0084327437    OMA13  NB7-C_18        C1        NbQ-B7
## 350   Otu6  WGA53 0.0082270671    WGA53 WGA_YE5-B        B1       NbP-YE5
## 43   Otu10  OMA18 0.0080213904    OMA18  RA5-B_18        A1        NbO-A5
## 200  Otu19  WGA18 0.0080213904    WGA18 WGA_RD4-A        A2        NbO-D4
## 250  Otu24  WGA54 0.0080213904    WGA54 WGA_YE5-B        B1       NbP-YE5
## 360   Otu6  WGA18 0.0080213904    WGA18 WGA_RD4-A        A2        NbO-D4
## 133  Otu12  WGA54 0.0078157137    WGA54 WGA_YE5-B        B1       NbP-YE5
## 124  Otu12  OMA13 0.0076100370    OMA13  NB7-C_18        C1        NbQ-B7
## 356   Otu6  WGA25 0.0076100370    WGA25 WGA_RD4-B        A2        NbO-D4
## 47   Otu10  WGA10 0.0074043603    WGA10 WGA_RA5-B        A1        NbO-A5
## 110  Otu12  OMA28 0.0074043603    OMA28  YE5-C_18        B1       NbP-YE5
## 116  Otu12  WGA46 0.0074043603    WGA46 WGA_YE5-A        B1       NbP-YE5
## 286   Otu4  OMA12 0.0074043603    OMA12  NB7-B_18        C1        NbQ-B7
## 112  Otu12   WGA3 0.0067873303     WGA3 WGA_RA5-A        A1        NbO-A5
## 274   Otu4  OMA11 0.0067873303    OMA11  NB7-A_18        C1        NbQ-B7
## 150  Otu14  OMA17 0.0065816536    OMA17  RA5-A_18        A1        NbO-A5
## 196  Otu19  WGA39 0.0065816536    WGA39 WGA_NB7-B        C1        NbQ-B7
## 55   Otu10  OMA19 0.0063759770    OMA19  RA5-C_18        A1        NbO-A5
## 45   Otu10  OMA17 0.0059646236    OMA17  RA5-A_18        A1        NbO-A5
## 156  Otu14  OMA26 0.0059646236    OMA26  YE5-A_18        B1       NbP-YE5
## 153  Otu14  OMA13 0.0057589469    OMA13  NB7-C_18        C1        NbQ-B7
## 160  Otu14  OMA20 0.0055532703    OMA20  RD4-A_18        A2        NbO-D4
## 254  Otu24   WGA3 0.0055532703     WGA3 WGA_RA5-A        A1        NbO-A5
## 287   Otu4  OMA28 0.0055532703    OMA28  YE5-C_18        B1       NbP-YE5
## 317   Otu5  WGA20 0.0055532703    WGA20 WGA_RD4-A        A2        NbO-D4
## 87   Otu11  OMA12 0.0049362402    OMA12  NB7-B_18        C1        NbQ-B7
## 234   Otu2  WGA41 0.0049362402    WGA41 WGA_NB7-B        C1        NbQ-B7
## 42   Otu10  OMA22 0.0047305636    OMA22  RD4-C_18        A2        NbO-D4
## 177  Otu19  OMA11 0.0047305636    OMA11  NB7-A_18        C1        NbQ-B7
## 182  Otu19  WGA20 0.0047305636    WGA20 WGA_RD4-A        A2        NbO-D4
## 375   Otu7  OMA11 0.0045248869    OMA11  NB7-A_18        C1        NbQ-B7
## 176  Otu19  OMA12 0.0043192102    OMA12  NB7-B_18        C1        NbQ-B7
## 404   Otu7  WGA18 0.0043192102    WGA18 WGA_RD4-A        A2        NbO-D4
## 102  Otu11  WGA55 0.0041135335    WGA55 WGA_YE5-B        B1       NbP-YE5
## 70   Otu11  OMA26 0.0037021802    OMA26  YE5-A_18        B1       NbP-YE5
## 158  Otu14  OMA21 0.0037021802    OMA21  RD4-B_18        A2        NbO-D4
## 383   Otu7  OMA27 0.0037021802    OMA27  YE5-B_18        B1       NbP-YE5
## 328   Otu5  WGA19 0.0034965035    WGA19 WGA_RD4-A        A2        NbO-D4
## 271  Otu24   WGA4 0.0032908268     WGA4 WGA_RA5-A        A1        NbO-A5
## 320   Otu5  WGA10 0.0032908268    WGA10 WGA_RA5-B        A1        NbO-A5
## 385   Otu7  OMA26 0.0032908268    OMA26  YE5-A_18        B1       NbP-YE5
## 239  Otu24  OMA22 0.0028794735    OMA22  RD4-C_18        A2        NbO-D4
## 108  Otu12  OMA20 0.0024681201    OMA20  RD4-A_18        A2        NbO-D4
## 258  Otu24  OMA20 0.0024681201    OMA20  RD4-A_18        A2        NbO-D4
## 381   Otu7  OMA12 0.0024681201    OMA12  NB7-B_18        C1        NbQ-B7
## 157  Otu14  OMA22 0.0022624434    OMA22  RD4-C_18        A2        NbO-D4
## 119  Otu12  OMA18 0.0018510901    OMA18  RA5-B_18        A1        NbO-A5
## 120  Otu12  OMA17 0.0016454134    OMA17  RA5-A_18        A1        NbO-A5
## 285   Otu4  OMA13 0.0014397367    OMA13  NB7-C_18        C1        NbQ-B7
## 380   Otu7  OMA13 0.0014397367    OMA13  NB7-C_18        C1        NbQ-B7
## 67   Otu10  WGA19 0.0012340601    WGA19 WGA_RD4-A        A2        NbO-D4
## 80   Otu11  WGA25 0.0012340601    WGA25 WGA_RD4-B        A2        NbO-D4
## 109  Otu12  OMA22 0.0012340601    OMA22  RD4-C_18        A2        NbO-D4
## 118  Otu12  OMA19 0.0012340601    OMA19  RA5-C_18        A1        NbO-A5
## 170  Otu14   WGA4 0.0012340601     WGA4 WGA_RA5-A        A1        NbO-A5
## 193  Otu19  WGA27 0.0012340601    WGA27 WGA_RD4-B        A2        NbO-D4
## 382   Otu7  OMA28 0.0012340601    OMA28  YE5-C_18        B1       NbP-YE5
## 255  Otu24   WGA5 0.0010283834     WGA5 WGA_RA5-A        A1        NbO-A5
## 105  Otu12  OMA21 0.0008227067    OMA21  RD4-B_18        A2        NbO-D4
## 389   Otu7  OMA20 0.0008227067    OMA20  RD4-A_18        A2        NbO-D4
## 57   Otu10  WGA13 0.0006170300    WGA13 WGA_RA5-B        A1        NbO-A5
## 304   Otu4  WGA39 0.0006170300    WGA39 WGA_NB7-B        C1        NbQ-B7
## 58   Otu10   WGA4 0.0004113534     WGA4 WGA_RA5-A        A1        NbO-A5
## 266  Otu24  WGA40 0.0004113534    WGA40 WGA_NB7-B        C1        NbQ-B7
## 272  Otu24  WGA39 0.0004113534    WGA39 WGA_NB7-B        C1        NbQ-B7
## 329   Otu5  WGA18 0.0004113534    WGA18 WGA_RD4-A        A2        NbO-D4
## 49   Otu10  WGA33 0.0002056767    WGA33 WGA_NB7-A        C1        NbQ-B7
## 64   Otu10   WGA5 0.0002056767     WGA5 WGA_RA5-A        A1        NbO-A5
## 65   Otu10  WGA41 0.0002056767    WGA41 WGA_NB7-B        C1        NbQ-B7
## 126  Otu12   WGA5 0.0002056767     WGA5 WGA_RA5-A        A1        NbO-A5
## 267  Otu24  WGA45 0.0002056767    WGA45 WGA_YE5-A        B1       NbP-YE5
## 348   Otu6  WGA13 0.0002056767    WGA13 WGA_RA5-B        A1        NbO-A5
## 39   Otu10  WGA25 0.0000000000    WGA25 WGA_RD4-B        A2        NbO-D4
## 44   Otu10  WGA20 0.0000000000    WGA20 WGA_RD4-A        A2        NbO-D4
## 50   Otu10   WGA3 0.0000000000     WGA3 WGA_RA5-A        A1        NbO-A5
## 54   Otu10  WGA45 0.0000000000    WGA45 WGA_YE5-A        B1       NbP-YE5
## 56   Otu10  WGA18 0.0000000000    WGA18 WGA_RD4-A        A2        NbO-D4
## 59   Otu10  WGA12 0.0000000000    WGA12 WGA_RA5-B        A1        NbO-A5
## 61   Otu10  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 66   Otu10  WGA40 0.0000000000    WGA40 WGA_NB7-B        C1        NbQ-B7
## 68   Otu10  WGA55 0.0000000000    WGA55 WGA_YE5-B        B1       NbP-YE5
## 89   Otu11  WGA33 0.0000000000    WGA33 WGA_NB7-A        C1        NbQ-B7
## 90   Otu11  WGA53 0.0000000000    WGA53 WGA_YE5-B        B1       NbP-YE5
## 91   Otu11  WGA45 0.0000000000    WGA45 WGA_YE5-A        B1       NbP-YE5
## 92   Otu11  WGA19 0.0000000000    WGA19 WGA_RD4-A        A2        NbO-D4
## 93   Otu11  WGA48 0.0000000000    WGA48 WGA_YE5-A        B1       NbP-YE5
## 94   Otu11  WGA46 0.0000000000    WGA46 WGA_YE5-A        B1       NbP-YE5
## 97   Otu11  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
## 100  Otu11  WGA54 0.0000000000    WGA54 WGA_YE5-B        B1       NbP-YE5
## 101  Otu11  WGA39 0.0000000000    WGA39 WGA_NB7-B        C1        NbQ-B7
## 106  Otu12  WGA25 0.0000000000    WGA25 WGA_RD4-B        A2        NbO-D4
## 107  Otu12  WGA20 0.0000000000    WGA20 WGA_RD4-A        A2        NbO-D4
## 114  Otu12  WGA27 0.0000000000    WGA27 WGA_RD4-B        A2        NbO-D4
## 121  Otu12  WGA12 0.0000000000    WGA12 WGA_RA5-B        A1        NbO-A5
## 122  Otu12  WGA10 0.0000000000    WGA10 WGA_RA5-B        A1        NbO-A5
## 127  Otu12  WGA19 0.0000000000    WGA19 WGA_RD4-A        A2        NbO-D4
## 128  Otu12  WGA18 0.0000000000    WGA18 WGA_RD4-A        A2        NbO-D4
## 129  Otu12  WGA13 0.0000000000    WGA13 WGA_RA5-B        A1        NbO-A5
## 130  Otu12   WGA4 0.0000000000     WGA4 WGA_RA5-A        A1        NbO-A5
## 142  Otu14  WGA20 0.0000000000    WGA20 WGA_RD4-A        A2        NbO-D4
## 143  Otu14  WGA46 0.0000000000    WGA46 WGA_YE5-A        B1       NbP-YE5
## 144  Otu14  WGA45 0.0000000000    WGA45 WGA_YE5-A        B1       NbP-YE5
## 147  Otu14  WGA18 0.0000000000    WGA18 WGA_RD4-A        A2        NbO-D4
## 149  Otu14  WGA12 0.0000000000    WGA12 WGA_RA5-B        A1        NbO-A5
## 151  Otu14  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 155  Otu14  WGA48 0.0000000000    WGA48 WGA_YE5-A        B1       NbP-YE5
## 163  Otu14  WGA39 0.0000000000    WGA39 WGA_NB7-B        C1        NbQ-B7
## 165  Otu14  WGA54 0.0000000000    WGA54 WGA_YE5-B        B1       NbP-YE5
## 166  Otu14  WGA53 0.0000000000    WGA53 WGA_YE5-B        B1       NbP-YE5
## 179  Otu19  WGA48 0.0000000000    WGA48 WGA_YE5-A        B1       NbP-YE5
## 183  Otu19  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
## 184  Otu19  WGA40 0.0000000000    WGA40 WGA_NB7-B        C1        NbQ-B7
## 188  Otu19  WGA53 0.0000000000    WGA53 WGA_YE5-B        B1       NbP-YE5
## 189  Otu19  WGA33 0.0000000000    WGA33 WGA_NB7-A        C1        NbQ-B7
## 195  Otu19  WGA55 0.0000000000    WGA55 WGA_YE5-B        B1       NbP-YE5
## 201  Otu19  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 246  Otu24  WGA27 0.0000000000    WGA27 WGA_RD4-B        A2        NbO-D4
## 247  Otu24  WGA25 0.0000000000    WGA25 WGA_RD4-B        A2        NbO-D4
## 251  Otu24  WGA53 0.0000000000    WGA53 WGA_YE5-B        B1       NbP-YE5
## 253  Otu24  WGA33 0.0000000000    WGA33 WGA_NB7-A        C1        NbQ-B7
## 257  Otu24  WGA46 0.0000000000    WGA46 WGA_YE5-A        B1       NbP-YE5
## 262  Otu24  WGA55 0.0000000000    WGA55 WGA_YE5-B        B1       NbP-YE5
## 268  Otu24  WGA19 0.0000000000    WGA19 WGA_RD4-A        A2        NbO-D4
## 269  Otu24  WGA18 0.0000000000    WGA18 WGA_RD4-A        A2        NbO-D4
## 270  Otu24  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
## 291   Otu4  WGA48 0.0000000000    WGA48 WGA_YE5-A        B1       NbP-YE5
## 300   Otu4  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 301   Otu4  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
## 321   Otu5  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 333   Otu5  WGA39 0.0000000000    WGA39 WGA_NB7-B        C1        NbQ-B7
## 340   Otu5  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
## 376   Otu7  WGA20 0.0000000000    WGA20 WGA_RD4-A        A2        NbO-D4
## 388   Otu7  WGA46 0.0000000000    WGA46 WGA_YE5-A        B1       NbP-YE5
## 395   Otu7  WGA34 0.0000000000    WGA34 WGA_NB7-A        C1        NbQ-B7
## 396   Otu7  WGA33 0.0000000000    WGA33 WGA_NB7-A        C1        NbQ-B7
## 397   Otu7  WGA53 0.0000000000    WGA53 WGA_YE5-B        B1       NbP-YE5
## 399   Otu7  WGA19 0.0000000000    WGA19 WGA_RD4-A        A2        NbO-D4
## 400   Otu7  WGA48 0.0000000000    WGA48 WGA_YE5-A        B1       NbP-YE5
## 401   Otu7  WGA40 0.0000000000    WGA40 WGA_NB7-B        C1        NbQ-B7
## 402   Otu7  WGA45 0.0000000000    WGA45 WGA_YE5-A        B1       NbP-YE5
## 403   Otu7  WGA39 0.0000000000    WGA39 WGA_NB7-B        C1        NbQ-B7
## 405   Otu7  WGA54 0.0000000000    WGA54 WGA_YE5-B        B1       NbP-YE5
## 407   Otu7  WGA55 0.0000000000    WGA55 WGA_YE5-B        B1       NbP-YE5
## 408   Otu7  WGA41 0.0000000000    WGA41 WGA_NB7-B        C1        NbQ-B7
##     Treatment Treatment3 Treatment4      Pop Pop2 Pop3 Sample_Site
## 296       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 295       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 281       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 111       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 18        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 224       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 231       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 222       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 7         NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 279       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 214       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 358       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 390       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 205       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 207       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 289       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 362       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 206       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 406       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 299       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 12        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 290       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 379       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 211       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 34        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 141       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 220       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 210       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 99        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 326       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 384       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 292       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 5         NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 13        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 256       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 305       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 181       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 277       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 28        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 1         RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 2         YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 30        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 32        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 168       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 21        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 235       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 4         YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 10        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 86        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 25        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 20        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 387       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 208       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 218       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 236       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 14        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 26        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 219       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 398       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 9         RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 288       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 24        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 276       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 228       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 223       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 22        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 3         RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 199       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 332       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 17        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 16        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 216       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 337       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 230       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 33        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 226       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 203       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 354       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 84        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 117       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 215       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 27        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 293       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 98        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 233       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 8         RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 131       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 306       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 186       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 29        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 213       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 134       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 95        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 238       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 322       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 393       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 391       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 343       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 227       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 297       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 15        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 148       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 237       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 370       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 217       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 394       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 209       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 313       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 232       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 174       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 187       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 6         RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 330       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 392       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 303       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 386       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 319       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 363       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 185       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 71        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 194       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 377       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 368       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 352       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 115       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 191       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 365       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 310       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 221       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 373       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 280       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 334       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 229       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 314       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 146       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 60        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 324       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 344       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 336       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 11        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 378       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 341       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 85        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 369       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 103       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 212       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 338       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 248       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 308       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 37        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 192       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 169       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 371       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 278       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 309       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 298       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 302       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 318       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 323       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 132       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 23        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 138       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 197       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 135       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 307       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 353       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 53        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 77        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 355       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 364       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 335       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 96        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 31        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 38        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 62        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 36        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 331       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 46        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 345       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 19        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 315       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 282       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 81        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 40        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 74        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 113       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 275       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 167       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 225       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 48        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 243       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 159       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 252       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 339       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 78        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 173       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 88        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 140       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 357       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 104       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 41        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 154       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 284       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 69        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 312       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 202       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 172       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 72        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 359       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 164       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 316       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 325       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 283       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 190       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 361       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 171       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 52        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 198       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 139       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 178       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 73        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 125       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 136       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 261       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 273       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 82        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 35        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 123       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 79        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 76        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 245       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 372       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 83        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 264       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 137       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 161       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 204       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 347       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 63        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 180       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 51        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 265       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 242       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 260       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 327       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 259       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 311       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 75        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 241       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 367       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 346       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 244       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 294       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 249       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 152       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 162       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 240       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 145       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 263       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 349       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 366       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 374       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 342       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 351       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 175       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 350       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 43        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 200       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 250       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 360       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 133       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 124       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 356       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 47        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 110       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 116       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 286       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 112       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 274       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 150       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 196       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 55        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 45        RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 156       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 153       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 160       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 254       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 287       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 317       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 87        NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 234       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 42        RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 177       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 182       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 375       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 176       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 404       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 102       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 70        YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 158       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 383       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 328       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 271       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 320       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 385       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 239       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 108       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 258       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 381       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 157       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 119       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 120       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 285       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 380       NB7        NB7        NB7      CG3  CG3 PopC         NbQ
## 67        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 80        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 109       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 118       RA5        RA5        RA5      CG1  CG1 PopA         NbO
## 170       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 193       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 382       YE5        YE5        YE5      CG2  CG2 PopB         NbP
## 255       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 105       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 389       RD4        RD4        RD4      CG1  CG1 PopA         NbO
## 57        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 304       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 58        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 266       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 272       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 329       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 49        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 64        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 65        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 126       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 267       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 348       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 39        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 44        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 50        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 54        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 56        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 59        RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 61        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 66        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 68        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 89        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 90        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 91        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 92        RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 93        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 94        YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 97        NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 100       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 101       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 106       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 107       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 114       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 121       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 122       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 127       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 128       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 129       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 130       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 142       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 143       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 144       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 147       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 149       RA5    RA5_WGA    RA5_WGA    March  CG1 PopA         NbO
## 151       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 155       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 163       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 165       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 166       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 179       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 183       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 184       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 188       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 189       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 195       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 201       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 246       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 247       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 251       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 253       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 257       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 262       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 268       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 269       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 270       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 291       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 300       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 301       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 321       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 333       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 340       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 376       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 388       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 395       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 396       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 397       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 399       RD4    RD4_WGA    RD4_WGA    March  CG1 PopA         NbO
## 400       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 401       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 402       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 403       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
## 405       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 407       YE5    YE5_WGA    YE5_WGA      May  CG2 PopB         NbP
## 408       NB7    NB7_WGA    NB7_WGA November  CG3 PopC         NbQ
##      treatment_26           treatment_26.1 extract         definition_26 temp
## 296  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 295  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 281     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 111  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 18   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 224  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 231     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 222     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 7    NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 279  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 214     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 358  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 390     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 205    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 207    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 289 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 362 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 206     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 406  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 299  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 12     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 290  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 379     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 211    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 34  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 141  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 220     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 210     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 99   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 326 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 384  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 292 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 5       NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 13      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 256 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 305 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 181 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 277     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 28     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 1    NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 2   NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 30   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 32      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 168  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 21      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 235  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 4      NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 10  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 86   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 25      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 20   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 387  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 208     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 218  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 236  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 14  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 26  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 219     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 398  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 9       NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 288  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 24   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 276  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 228 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 223  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 22  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 3    NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 199  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 332  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 17      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 16   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 216  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 337  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 230  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 33      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 226 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 203  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 354 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 84   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 117 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 215  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 27   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 293  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 98   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 233  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 8    NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 131  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 306  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 186     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 29      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 213  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 134 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 95   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 238  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 322 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 393  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 391     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 343 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 227  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 297 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 15   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 148  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 237 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 370 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 217 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 394  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 209     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 313  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 232  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 174     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 187 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 6    NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 330  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 392  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 303  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 386  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 319     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 363  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 185  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 71      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 194  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 377     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 368    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 352    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 115 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 191  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 365     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 310    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 221  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 373     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 280 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 334 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 229 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 314    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 146     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 60   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 324  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 344  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 336 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 11   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 378     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 341     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 85   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 369    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 103    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 212  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 338 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 248  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 308    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 37     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 192  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 169  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 371  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 278  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 309     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 298  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 302  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 318     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 323  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 132 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 23   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 138  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 197     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 135  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 307     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 353     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 53  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 77   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 355  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 364  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 335  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 96   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 31   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 38      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 62  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 36      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 331     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 46      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 345     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 19   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 315     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 282     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 81   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 40     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 74      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 113     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 275    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 167  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 225 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 48     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 243     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 159  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 252    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 339  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 78   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 173     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 88     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 140    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 357     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 104    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 41      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 154     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 284     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 69      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 312     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 202  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 172     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 72      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 359  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 164 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 316     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 325  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 283     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 190    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 361     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 171    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 52  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 198     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 139  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 178    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 73      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 125     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 136  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 261  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 273    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 82      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 35      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 123  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 79   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 76     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 245    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 372  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 83      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 264  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 137  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 161     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 204  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 347  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 63  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 180 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 51   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 265  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 242     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 260     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 327 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 259     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 311     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 75      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 241     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 367  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 346  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 244    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 294     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 249     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 152  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 162     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 240     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 145    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 263  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 349     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 366  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 374     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 342 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 351  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 175     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 350 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 43      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 200  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 250 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 360  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 133 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 124     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 356  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 47   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 110    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 116 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 286     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 112  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 274     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 150     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 196  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 55      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 45      NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 156    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 153     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 160     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 254  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 287    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 317  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 87      NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 234  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 42      NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 177     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 182  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 375     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 176     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 404  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 102 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 70     NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 158     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 383    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 328  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 271  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 320  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 385    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 239     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 108     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 258     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 381     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 157     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 119     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 120     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 285     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 380     NbQ-B7 BA             C1 Coculture  Filter             Coculture   18
## 67   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 80   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 109     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 118     NbO-A5 BA             A1 Coculture  Filter             Coculture   18
## 170  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 193  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 382    NbP-YE5 BA             B1 Coculture  Filter             Coculture   18
## 255  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 105     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 389     NbO-D4 BA             A2 Coculture  Filter             Coculture   18
## 57   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 304  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 58   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 266  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 272  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 329  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 49   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 64   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 65   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 126  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 267 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 348  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 39   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 44   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 50   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 54  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 56   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 59   NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 61   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 66   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 68  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 89   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 90  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 91  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 92   NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 93  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 94  NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 97   NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 100 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 101  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 106  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 107  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 114  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 121  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 122  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 127  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 128  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 129  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 130  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 142  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 143 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 144 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 147  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 149  NbO-A5 Micro A1 Algal-cell associated     WGA Algal-cell associated   18
## 151  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 155 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 163  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 165 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 166 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 179 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 183  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 184  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 188 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 189  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 195 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 201  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 246  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 247  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 251 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 253  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 257 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 262 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 268  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 269  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 270  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 291 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 300  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 301  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 321  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 333  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 340  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 376  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 388 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 395  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 396  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 397 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 399  NbO-D4 Micro A2 Algal-cell associated     WGA Algal-cell associated   18
## 400 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 401  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 402 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 403  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
## 405 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 407 NbP-YE5 Micro B1 Algal-cell associated     WGA Algal-cell associated   18
## 408  NbQ-B7 Micro C1 Algal-cell associated     WGA Algal-cell associated   18
##             col    col2  col_26 pch culture definition keep_2026 diatom_control
## 296        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 295        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 281        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 111 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 18  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 224 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 231        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 222 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 7   forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 279        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 214 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 358 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 390        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 205   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 207   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 289   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 362   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 206 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 406  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 299  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 12    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 290        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 379        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 211   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 34    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 141        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 220        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 210        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 99  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 326   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 384  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 292   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 5   forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 13  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 256   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 305   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 181   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 277        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 28    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 1    steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 2     firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 30  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 32         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 168 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 21         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 235  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 4     firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 10    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 86   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 25  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 20   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 387        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 208        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 218        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 236 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 14    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 26    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 219        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 398  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 9          navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 288 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 24   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 276  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 228   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 223  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 22    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 3    steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 199        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 332  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 17         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 16         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 216        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 337  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 230 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 33         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 226   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 203  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 354   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 84         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 117   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 215        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 27         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 293  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 98  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 233        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 8    steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 131 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 306  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 186        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 29         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 213  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 134   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 95   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 238  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 322   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 393  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 391        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 343   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 227  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 297   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 15         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 148  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 237   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 370   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 217   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 394  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 209        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 313 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 232        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 174        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 187   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 6    steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 330  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 392  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 303  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 386        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 319 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 363 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 185  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 71         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 194        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 377        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 368   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 352   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 115   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 191  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 365 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 310   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 221  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 373        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 280   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 334   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 229   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 314   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 146 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 60  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 324        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 344        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 336   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 11  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 378        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 341 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 85   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 369   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 103   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 212 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 338   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 248        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 308   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 37    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 192  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 169 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 371 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 278        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 309 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 298  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 302 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 318        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 323  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 132   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 23  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 138        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 197        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 135 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 307 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 353 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 53    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 77   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 355        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 364  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 335 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 96   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 31         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 38  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 62    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 36  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 331        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 46  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 345        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 19         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 315        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 282        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 81         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 40    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 74         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 113 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 275   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 167  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 225   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 48    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 243 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 159        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 252   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 339  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 78   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 173        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 88    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 140   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 357        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 104   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 41         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 154 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 284        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 69         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 312        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 202  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 172        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 72         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 359        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 164   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 316        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 325        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 283        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 190   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 361        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 171   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 52    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 198        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 139 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 178   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 73  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 125 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 136 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 261  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 273   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 82         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 35         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 123 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 79         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 76    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 245   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 372  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 83         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 264  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 137  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 161        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 204  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 347  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 63    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 180   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 51         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 265 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 242 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 260        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 327   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 259        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 311        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 75  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 241 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 367  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 346 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 244   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 294        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 249        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 152  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 162        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 240        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 145   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 263  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 349        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 366  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 374        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 342   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 351 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 175 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 350   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 43         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 200        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 250   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 360        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 133   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 124 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 356        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 47   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 110   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 116   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 286 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 112  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 274 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 150        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 196 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 55         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 45         navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 156   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 153 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 160        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 254  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 287   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 317        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 87  forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 234 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 42         navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 177 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 182        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 375 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 176 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 404        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 102   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 70    firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 158        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 383   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 328        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 271  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 320  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 385   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 239        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 108        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 258        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 381 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 157        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 119        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 120        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 285 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 380 forestgreen #7768AE #3BB273  21  Single assemblage       yes         diatom
## 67         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 80         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 109        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 118        navy #4D9DE0 #7768AE  21  Single assemblage       yes         diatom
## 170  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 193        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 382   firebrick #3BB273 #E15554  21  Single assemblage       yes         diatom
## 255  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 105        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 389        navy #E15554 #4D9DE0  21  Single assemblage       yes         diatom
## 57   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 304 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 58   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 266 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 272 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 329        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 49  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 64   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 65  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 126  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 267   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 348  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 39         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 44         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 50   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 54    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 56         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 59   steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 61  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 66  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 68    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 89  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 90    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 91    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 92         navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 93    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 94    firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 97  forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 100   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 101 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 106        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 107        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 114        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 121  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 122  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 127        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 128        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 129  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 130  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 142        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 143   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 144   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 147        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 149  steelblue4 #4D9DE0 #7768AE  22  Single microbiome       yes         diatom
## 151 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 155   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 163 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 165   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 166   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 179   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 183 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 184 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 188   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 189 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 195   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 201 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 246        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 247        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 251   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 253 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 257   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 262   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 268        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 269        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 270 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 291   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 300 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 301 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 321 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 333 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 340 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 376        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 388   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 395 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 396 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 397   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 399        navy #E15554 #4D9DE0  22  Single microbiome       yes         diatom
## 400   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 401 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 402   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 403 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
## 405   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 407   firebrick #3BB273 #E15554  22  Single microbiome       yes         diatom
## 408 forestgreen #7768AE #3BB273  22  Single microbiome       yes         diatom
##     innoc  Kingdom          phyla               class            family
## 296  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 295  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 281  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 111  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 18   <NA>     <NA>           <NA>                <NA>              <NA>
## 224  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 231  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 222  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 7    <NA>     <NA>           <NA>                <NA>              <NA>
## 279  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 214  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 358  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 390  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 205  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 207  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 289  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 362  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 206  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 406  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 299  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 12   <NA>     <NA>           <NA>                <NA>              <NA>
## 290  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 379  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 211  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 34   <NA>     <NA>           <NA>                <NA>              <NA>
## 141  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 220  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 210  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 99   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 326  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 384  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 292  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 5    <NA>     <NA>           <NA>                <NA>              <NA>
## 13   <NA>     <NA>           <NA>                <NA>              <NA>
## 256  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 305  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 181  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 277  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 28   <NA>     <NA>           <NA>                <NA>              <NA>
## 1    <NA>     <NA>           <NA>                <NA>              <NA>
## 2    <NA>     <NA>           <NA>                <NA>              <NA>
## 30   <NA>     <NA>           <NA>                <NA>              <NA>
## 32   <NA>     <NA>           <NA>                <NA>              <NA>
## 168  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 21   <NA>     <NA>           <NA>                <NA>              <NA>
## 235  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 4    <NA>     <NA>           <NA>                <NA>              <NA>
## 10   <NA>     <NA>           <NA>                <NA>              <NA>
## 86   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 25   <NA>     <NA>           <NA>                <NA>              <NA>
## 20   <NA>     <NA>           <NA>                <NA>              <NA>
## 387  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 208  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 218  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 236  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 14   <NA>     <NA>           <NA>                <NA>              <NA>
## 26   <NA>     <NA>           <NA>                <NA>              <NA>
## 219  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 398  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 9    <NA>     <NA>           <NA>                <NA>              <NA>
## 288  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 24   <NA>     <NA>           <NA>                <NA>              <NA>
## 276  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 228  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 223  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 22   <NA>     <NA>           <NA>                <NA>              <NA>
## 3    <NA>     <NA>           <NA>                <NA>              <NA>
## 199  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 332  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 17   <NA>     <NA>           <NA>                <NA>              <NA>
## 16   <NA>     <NA>           <NA>                <NA>              <NA>
## 216  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 337  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 230  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 33   <NA>     <NA>           <NA>                <NA>              <NA>
## 226  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 203  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 354  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 84   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 117  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 215  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 27   <NA>     <NA>           <NA>                <NA>              <NA>
## 293  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 98   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 233  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 8    <NA>     <NA>           <NA>                <NA>              <NA>
## 131  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 306  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 186  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 29   <NA>     <NA>           <NA>                <NA>              <NA>
## 213  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 134  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 95   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 238  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 322  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 393  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 391  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 343  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 227  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 297  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 15   <NA>     <NA>           <NA>                <NA>              <NA>
## 148  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 237  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 370  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 217  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 394  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 209  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 313  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 232  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 174  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 187  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 6    <NA>     <NA>           <NA>                <NA>              <NA>
## 330  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 392  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 303  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 386  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 319  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 363  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 185  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 71   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 194  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 377  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 368  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 352  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 115  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 191  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 365  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 310  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 221  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 373  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 280  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 334  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 229  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 314  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 146  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 60   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 324  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 344  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 336  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 11   <NA>     <NA>           <NA>                <NA>              <NA>
## 378  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 341  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 85   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 369  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 103  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 212  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 338  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 248  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 308  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 37   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 192  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 169  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 371  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 278  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 309  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 298  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 302  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 318  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 323  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 132  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 23   <NA>     <NA>           <NA>                <NA>              <NA>
## 138  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 197  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 135  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 307  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 353  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 53   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 77   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 355  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 364  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 335  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 96   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 31   <NA>     <NA>           <NA>                <NA>              <NA>
## 38   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 62   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 36   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 331  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 46   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 345  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 19   <NA>     <NA>           <NA>                <NA>              <NA>
## 315  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 282  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 81   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 40   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 74   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 113  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 275  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 167  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 225  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 48   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 243  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 159  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 252  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 339  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 78   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 173  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 88   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 140  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 357  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 104  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 41   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 154  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 284  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 69   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 312  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 202  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 172  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 72   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 359  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 164  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 316  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 325  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 283  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 190  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 361  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 171  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 52   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 198  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 139  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 178  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 73   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 125  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 136  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 261  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 273  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 82   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 35   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 123  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 79   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 76   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 245  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 372  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 83   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 264  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 137  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 161  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 204  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 347  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 63   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 180  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 51   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 265  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 242  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 260  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 327  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 259  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 311  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 75   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 241  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 367  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 346  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 244  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 294  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 249  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 152  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 162  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 240  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 145  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 263  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 349  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 366  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 374  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 342  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 351  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 175  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 350  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 43   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 200  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 250  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 360  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 133  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 124  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 356  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 47   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 110  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 116  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 286  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 112  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 274  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 150  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 196  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 55   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 45   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 156  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 153  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 160  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 254  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 287  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 317  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 87   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 234  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 42   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 177  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 182  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 375  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 176  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 404  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 102  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 70   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 158  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 383  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 328  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 271  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 320  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 385  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 239  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 108  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 258  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 381  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 157  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 119  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 120  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 285  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 380  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 67   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 80   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 109  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 118  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 170  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 193  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 382  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 255  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 105  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 389  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 57   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 304  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 58   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 266  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 272  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 329  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 49   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 64   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 65   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 126  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 267  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 348  <NA> Bacteria Proteobacteria Alphaproteobacteria             SAR11
## 39   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 44   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 50   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 54   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 56   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 59   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 61   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 66   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 68   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 89   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 90   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 91   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 92   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 93   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 94   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 97   <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 100  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 101  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 106  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 107  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 114  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 121  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 122  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 127  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 128  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 129  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 130  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 142  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 143  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 144  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 147  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 149  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 151  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 155  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 163  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 165  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 166  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 179  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 183  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 184  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 188  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 189  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 195  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 201  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 246  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 247  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 251  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 253  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 257  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 262  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 268  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 269  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 270  <NA> Bacteria  Bacteroidetes      Flavobacteriia  Flavobacteriales
## 291  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 300  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 301  <NA> Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 321  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 333  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 340  <NA> Bacteria Proteobacteria Gammaproteobacteria   Alteromonadales
## 376  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 388  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 395  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 396  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 397  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 399  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 400  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 401  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 402  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 403  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 405  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 407  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
## 408  <NA> Bacteria Proteobacteria Gammaproteobacteria       Vibrionales
##                  genus                 species
## 296 Oceanospirillaceae             Marinomonas
## 295 Oceanospirillaceae             Marinomonas
## 281 Oceanospirillaceae             Marinomonas
## 111 Oceanospirillaceae          Neptuniibacter
## 18                <NA>                    <NA>
## 224  Flavobacteriaceae           Tenacibaculum
## 231  Flavobacteriaceae           Tenacibaculum
## 222  Flavobacteriaceae           Tenacibaculum
## 7                 <NA>                    <NA>
## 279 Oceanospirillaceae             Marinomonas
## 214  Flavobacteriaceae           Tenacibaculum
## 358              SAR11 Candidatus Pelagibacter
## 390       Vibrionaceae                  Vibrio
## 205  Flavobacteriaceae           Tenacibaculum
## 207  Flavobacteriaceae           Tenacibaculum
## 289 Oceanospirillaceae             Marinomonas
## 362              SAR11 Candidatus Pelagibacter
## 206  Flavobacteriaceae           Tenacibaculum
## 406       Vibrionaceae                  Vibrio
## 299 Oceanospirillaceae             Marinomonas
## 12                <NA>                    <NA>
## 290 Oceanospirillaceae             Marinomonas
## 379       Vibrionaceae                  Vibrio
## 211  Flavobacteriaceae           Tenacibaculum
## 34                <NA>                    <NA>
## 141 Oceanospirillaceae         Marinobacterium
## 220  Flavobacteriaceae           Tenacibaculum
## 210  Flavobacteriaceae           Tenacibaculum
## 99      Cryomorphaceae             Owenweeksia
## 326   Alteromonadaceae         Aestuariibacter
## 384       Vibrionaceae                  Vibrio
## 292 Oceanospirillaceae             Marinomonas
## 5                 <NA>                    <NA>
## 13                <NA>                    <NA>
## 256  Flavobacteriaceae              Aquibacter
## 305 Oceanospirillaceae             Marinomonas
## 181     Cryomorphaceae                Wandonia
## 277 Oceanospirillaceae             Marinomonas
## 28                <NA>                    <NA>
## 1                 <NA>                    <NA>
## 2                 <NA>                    <NA>
## 30                <NA>                    <NA>
## 32                <NA>                    <NA>
## 168 Oceanospirillaceae         Marinobacterium
## 21                <NA>                    <NA>
## 235  Flavobacteriaceae           Tenacibaculum
## 4                 <NA>                    <NA>
## 10                <NA>                    <NA>
## 86      Cryomorphaceae             Owenweeksia
## 25                <NA>                    <NA>
## 20                <NA>                    <NA>
## 387       Vibrionaceae                  Vibrio
## 208  Flavobacteriaceae           Tenacibaculum
## 218  Flavobacteriaceae           Tenacibaculum
## 236  Flavobacteriaceae           Tenacibaculum
## 14                <NA>                    <NA>
## 26                <NA>                    <NA>
## 219  Flavobacteriaceae           Tenacibaculum
## 398       Vibrionaceae                  Vibrio
## 9                 <NA>                    <NA>
## 288 Oceanospirillaceae             Marinomonas
## 24                <NA>                    <NA>
## 276 Oceanospirillaceae             Marinomonas
## 228  Flavobacteriaceae           Tenacibaculum
## 223  Flavobacteriaceae           Tenacibaculum
## 22                <NA>                    <NA>
## 3                 <NA>                    <NA>
## 199     Cryomorphaceae                Wandonia
## 332   Alteromonadaceae         Aestuariibacter
## 17                <NA>                    <NA>
## 16                <NA>                    <NA>
## 216  Flavobacteriaceae           Tenacibaculum
## 337   Alteromonadaceae         Aestuariibacter
## 230  Flavobacteriaceae           Tenacibaculum
## 33                <NA>                    <NA>
## 226  Flavobacteriaceae           Tenacibaculum
## 203     Cryomorphaceae                Wandonia
## 354              SAR11 Candidatus Pelagibacter
## 84      Cryomorphaceae             Owenweeksia
## 117 Oceanospirillaceae          Neptuniibacter
## 215  Flavobacteriaceae           Tenacibaculum
## 27                <NA>                    <NA>
## 293 Oceanospirillaceae             Marinomonas
## 98      Cryomorphaceae             Owenweeksia
## 233  Flavobacteriaceae           Tenacibaculum
## 8                 <NA>                    <NA>
## 131 Oceanospirillaceae          Neptuniibacter
## 306 Oceanospirillaceae             Marinomonas
## 186     Cryomorphaceae                Wandonia
## 29                <NA>                    <NA>
## 213  Flavobacteriaceae           Tenacibaculum
## 134 Oceanospirillaceae          Neptuniibacter
## 95      Cryomorphaceae             Owenweeksia
## 238  Flavobacteriaceae           Tenacibaculum
## 322   Alteromonadaceae         Aestuariibacter
## 393       Vibrionaceae                  Vibrio
## 391       Vibrionaceae                  Vibrio
## 343              SAR11 Candidatus Pelagibacter
## 227  Flavobacteriaceae           Tenacibaculum
## 297 Oceanospirillaceae             Marinomonas
## 15                <NA>                    <NA>
## 148 Oceanospirillaceae         Marinobacterium
## 237  Flavobacteriaceae           Tenacibaculum
## 370              SAR11 Candidatus Pelagibacter
## 217  Flavobacteriaceae           Tenacibaculum
## 394       Vibrionaceae                  Vibrio
## 209  Flavobacteriaceae           Tenacibaculum
## 313   Alteromonadaceae         Aestuariibacter
## 232  Flavobacteriaceae           Tenacibaculum
## 174     Cryomorphaceae                Wandonia
## 187     Cryomorphaceae                Wandonia
## 6                 <NA>                    <NA>
## 330   Alteromonadaceae         Aestuariibacter
## 392       Vibrionaceae                  Vibrio
## 303 Oceanospirillaceae             Marinomonas
## 386       Vibrionaceae                  Vibrio
## 319   Alteromonadaceae         Aestuariibacter
## 363              SAR11 Candidatus Pelagibacter
## 185     Cryomorphaceae                Wandonia
## 71      Cryomorphaceae             Owenweeksia
## 194     Cryomorphaceae                Wandonia
## 377       Vibrionaceae                  Vibrio
## 368              SAR11 Candidatus Pelagibacter
## 352              SAR11 Candidatus Pelagibacter
## 115 Oceanospirillaceae          Neptuniibacter
## 191     Cryomorphaceae                Wandonia
## 365              SAR11 Candidatus Pelagibacter
## 310   Alteromonadaceae         Aestuariibacter
## 221  Flavobacteriaceae           Tenacibaculum
## 373              SAR11 Candidatus Pelagibacter
## 280 Oceanospirillaceae             Marinomonas
## 334   Alteromonadaceae         Aestuariibacter
## 229  Flavobacteriaceae           Tenacibaculum
## 314   Alteromonadaceae         Aestuariibacter
## 146 Oceanospirillaceae         Marinobacterium
## 60   Flavobacteriaceae            Cellulophaga
## 324   Alteromonadaceae         Aestuariibacter
## 344              SAR11 Candidatus Pelagibacter
## 336   Alteromonadaceae         Aestuariibacter
## 11                <NA>                    <NA>
## 378       Vibrionaceae                  Vibrio
## 341              SAR11 Candidatus Pelagibacter
## 85      Cryomorphaceae             Owenweeksia
## 369              SAR11 Candidatus Pelagibacter
## 103 Oceanospirillaceae          Neptuniibacter
## 212  Flavobacteriaceae           Tenacibaculum
## 338   Alteromonadaceae         Aestuariibacter
## 248  Flavobacteriaceae              Aquibacter
## 308   Alteromonadaceae         Aestuariibacter
## 37   Flavobacteriaceae            Cellulophaga
## 192     Cryomorphaceae                Wandonia
## 169 Oceanospirillaceae         Marinobacterium
## 371              SAR11 Candidatus Pelagibacter
## 278 Oceanospirillaceae             Marinomonas
## 309   Alteromonadaceae         Aestuariibacter
## 298 Oceanospirillaceae             Marinomonas
## 302 Oceanospirillaceae             Marinomonas
## 318   Alteromonadaceae         Aestuariibacter
## 323   Alteromonadaceae         Aestuariibacter
## 132 Oceanospirillaceae          Neptuniibacter
## 23                <NA>                    <NA>
## 138 Oceanospirillaceae         Marinobacterium
## 197     Cryomorphaceae                Wandonia
## 135 Oceanospirillaceae          Neptuniibacter
## 307   Alteromonadaceae         Aestuariibacter
## 353              SAR11 Candidatus Pelagibacter
## 53   Flavobacteriaceae            Cellulophaga
## 77      Cryomorphaceae             Owenweeksia
## 355              SAR11 Candidatus Pelagibacter
## 364              SAR11 Candidatus Pelagibacter
## 335   Alteromonadaceae         Aestuariibacter
## 96      Cryomorphaceae             Owenweeksia
## 31                <NA>                    <NA>
## 38   Flavobacteriaceae            Cellulophaga
## 62   Flavobacteriaceae            Cellulophaga
## 36   Flavobacteriaceae            Cellulophaga
## 331   Alteromonadaceae         Aestuariibacter
## 46   Flavobacteriaceae            Cellulophaga
## 345              SAR11 Candidatus Pelagibacter
## 19                <NA>                    <NA>
## 315   Alteromonadaceae         Aestuariibacter
## 282 Oceanospirillaceae             Marinomonas
## 81      Cryomorphaceae             Owenweeksia
## 40   Flavobacteriaceae            Cellulophaga
## 74      Cryomorphaceae             Owenweeksia
## 113 Oceanospirillaceae          Neptuniibacter
## 275 Oceanospirillaceae             Marinomonas
## 167 Oceanospirillaceae         Marinobacterium
## 225  Flavobacteriaceae           Tenacibaculum
## 48   Flavobacteriaceae            Cellulophaga
## 243  Flavobacteriaceae              Aquibacter
## 159 Oceanospirillaceae         Marinobacterium
## 252  Flavobacteriaceae              Aquibacter
## 339   Alteromonadaceae         Aestuariibacter
## 78      Cryomorphaceae             Owenweeksia
## 173     Cryomorphaceae                Wandonia
## 88      Cryomorphaceae             Owenweeksia
## 140 Oceanospirillaceae         Marinobacterium
## 357              SAR11 Candidatus Pelagibacter
## 104 Oceanospirillaceae          Neptuniibacter
## 41   Flavobacteriaceae            Cellulophaga
## 154 Oceanospirillaceae         Marinobacterium
## 284 Oceanospirillaceae             Marinomonas
## 69      Cryomorphaceae             Owenweeksia
## 312   Alteromonadaceae         Aestuariibacter
## 202     Cryomorphaceae                Wandonia
## 172     Cryomorphaceae                Wandonia
## 72      Cryomorphaceae             Owenweeksia
## 359              SAR11 Candidatus Pelagibacter
## 164 Oceanospirillaceae         Marinobacterium
## 316   Alteromonadaceae         Aestuariibacter
## 325   Alteromonadaceae         Aestuariibacter
## 283 Oceanospirillaceae             Marinomonas
## 190     Cryomorphaceae                Wandonia
## 361              SAR11 Candidatus Pelagibacter
## 171     Cryomorphaceae                Wandonia
## 52   Flavobacteriaceae            Cellulophaga
## 198     Cryomorphaceae                Wandonia
## 139 Oceanospirillaceae         Marinobacterium
## 178     Cryomorphaceae                Wandonia
## 73      Cryomorphaceae             Owenweeksia
## 125 Oceanospirillaceae          Neptuniibacter
## 136 Oceanospirillaceae          Neptuniibacter
## 261  Flavobacteriaceae              Aquibacter
## 273 Oceanospirillaceae             Marinomonas
## 82      Cryomorphaceae             Owenweeksia
## 35   Flavobacteriaceae            Cellulophaga
## 123 Oceanospirillaceae          Neptuniibacter
## 79      Cryomorphaceae             Owenweeksia
## 76      Cryomorphaceae             Owenweeksia
## 245  Flavobacteriaceae              Aquibacter
## 372              SAR11 Candidatus Pelagibacter
## 83      Cryomorphaceae             Owenweeksia
## 264  Flavobacteriaceae              Aquibacter
## 137 Oceanospirillaceae         Marinobacterium
## 161 Oceanospirillaceae         Marinobacterium
## 204     Cryomorphaceae                Wandonia
## 347              SAR11 Candidatus Pelagibacter
## 63   Flavobacteriaceae            Cellulophaga
## 180     Cryomorphaceae                Wandonia
## 51   Flavobacteriaceae            Cellulophaga
## 265  Flavobacteriaceae              Aquibacter
## 242  Flavobacteriaceae              Aquibacter
## 260  Flavobacteriaceae              Aquibacter
## 327   Alteromonadaceae         Aestuariibacter
## 259  Flavobacteriaceae              Aquibacter
## 311   Alteromonadaceae         Aestuariibacter
## 75      Cryomorphaceae             Owenweeksia
## 241  Flavobacteriaceae              Aquibacter
## 367              SAR11 Candidatus Pelagibacter
## 346              SAR11 Candidatus Pelagibacter
## 244  Flavobacteriaceae              Aquibacter
## 294 Oceanospirillaceae             Marinomonas
## 249  Flavobacteriaceae              Aquibacter
## 152 Oceanospirillaceae         Marinobacterium
## 162 Oceanospirillaceae         Marinobacterium
## 240  Flavobacteriaceae              Aquibacter
## 145 Oceanospirillaceae         Marinobacterium
## 263  Flavobacteriaceae              Aquibacter
## 349              SAR11 Candidatus Pelagibacter
## 366              SAR11 Candidatus Pelagibacter
## 374              SAR11 Candidatus Pelagibacter
## 342              SAR11 Candidatus Pelagibacter
## 351              SAR11 Candidatus Pelagibacter
## 175     Cryomorphaceae                Wandonia
## 350              SAR11 Candidatus Pelagibacter
## 43   Flavobacteriaceae            Cellulophaga
## 200     Cryomorphaceae                Wandonia
## 250  Flavobacteriaceae              Aquibacter
## 360              SAR11 Candidatus Pelagibacter
## 133 Oceanospirillaceae          Neptuniibacter
## 124 Oceanospirillaceae          Neptuniibacter
## 356              SAR11 Candidatus Pelagibacter
## 47   Flavobacteriaceae            Cellulophaga
## 110 Oceanospirillaceae          Neptuniibacter
## 116 Oceanospirillaceae          Neptuniibacter
## 286 Oceanospirillaceae             Marinomonas
## 112 Oceanospirillaceae          Neptuniibacter
## 274 Oceanospirillaceae             Marinomonas
## 150 Oceanospirillaceae         Marinobacterium
## 196     Cryomorphaceae                Wandonia
## 55   Flavobacteriaceae            Cellulophaga
## 45   Flavobacteriaceae            Cellulophaga
## 156 Oceanospirillaceae         Marinobacterium
## 153 Oceanospirillaceae         Marinobacterium
## 160 Oceanospirillaceae         Marinobacterium
## 254  Flavobacteriaceae              Aquibacter
## 287 Oceanospirillaceae             Marinomonas
## 317   Alteromonadaceae         Aestuariibacter
## 87      Cryomorphaceae             Owenweeksia
## 234  Flavobacteriaceae           Tenacibaculum
## 42   Flavobacteriaceae            Cellulophaga
## 177     Cryomorphaceae                Wandonia
## 182     Cryomorphaceae                Wandonia
## 375       Vibrionaceae                  Vibrio
## 176     Cryomorphaceae                Wandonia
## 404       Vibrionaceae                  Vibrio
## 102     Cryomorphaceae             Owenweeksia
## 70      Cryomorphaceae             Owenweeksia
## 158 Oceanospirillaceae         Marinobacterium
## 383       Vibrionaceae                  Vibrio
## 328   Alteromonadaceae         Aestuariibacter
## 271  Flavobacteriaceae              Aquibacter
## 320   Alteromonadaceae         Aestuariibacter
## 385       Vibrionaceae                  Vibrio
## 239  Flavobacteriaceae              Aquibacter
## 108 Oceanospirillaceae          Neptuniibacter
## 258  Flavobacteriaceae              Aquibacter
## 381       Vibrionaceae                  Vibrio
## 157 Oceanospirillaceae         Marinobacterium
## 119 Oceanospirillaceae          Neptuniibacter
## 120 Oceanospirillaceae          Neptuniibacter
## 285 Oceanospirillaceae             Marinomonas
## 380       Vibrionaceae                  Vibrio
## 67   Flavobacteriaceae            Cellulophaga
## 80      Cryomorphaceae             Owenweeksia
## 109 Oceanospirillaceae          Neptuniibacter
## 118 Oceanospirillaceae          Neptuniibacter
## 170 Oceanospirillaceae         Marinobacterium
## 193     Cryomorphaceae                Wandonia
## 382       Vibrionaceae                  Vibrio
## 255  Flavobacteriaceae              Aquibacter
## 105 Oceanospirillaceae          Neptuniibacter
## 389       Vibrionaceae                  Vibrio
## 57   Flavobacteriaceae            Cellulophaga
## 304 Oceanospirillaceae             Marinomonas
## 58   Flavobacteriaceae            Cellulophaga
## 266  Flavobacteriaceae              Aquibacter
## 272  Flavobacteriaceae              Aquibacter
## 329   Alteromonadaceae         Aestuariibacter
## 49   Flavobacteriaceae            Cellulophaga
## 64   Flavobacteriaceae            Cellulophaga
## 65   Flavobacteriaceae            Cellulophaga
## 126 Oceanospirillaceae          Neptuniibacter
## 267  Flavobacteriaceae              Aquibacter
## 348              SAR11 Candidatus Pelagibacter
## 39   Flavobacteriaceae            Cellulophaga
## 44   Flavobacteriaceae            Cellulophaga
## 50   Flavobacteriaceae            Cellulophaga
## 54   Flavobacteriaceae            Cellulophaga
## 56   Flavobacteriaceae            Cellulophaga
## 59   Flavobacteriaceae            Cellulophaga
## 61   Flavobacteriaceae            Cellulophaga
## 66   Flavobacteriaceae            Cellulophaga
## 68   Flavobacteriaceae            Cellulophaga
## 89      Cryomorphaceae             Owenweeksia
## 90      Cryomorphaceae             Owenweeksia
## 91      Cryomorphaceae             Owenweeksia
## 92      Cryomorphaceae             Owenweeksia
## 93      Cryomorphaceae             Owenweeksia
## 94      Cryomorphaceae             Owenweeksia
## 97      Cryomorphaceae             Owenweeksia
## 100     Cryomorphaceae             Owenweeksia
## 101     Cryomorphaceae             Owenweeksia
## 106 Oceanospirillaceae          Neptuniibacter
## 107 Oceanospirillaceae          Neptuniibacter
## 114 Oceanospirillaceae          Neptuniibacter
## 121 Oceanospirillaceae          Neptuniibacter
## 122 Oceanospirillaceae          Neptuniibacter
## 127 Oceanospirillaceae          Neptuniibacter
## 128 Oceanospirillaceae          Neptuniibacter
## 129 Oceanospirillaceae          Neptuniibacter
## 130 Oceanospirillaceae          Neptuniibacter
## 142 Oceanospirillaceae         Marinobacterium
## 143 Oceanospirillaceae         Marinobacterium
## 144 Oceanospirillaceae         Marinobacterium
## 147 Oceanospirillaceae         Marinobacterium
## 149 Oceanospirillaceae         Marinobacterium
## 151 Oceanospirillaceae         Marinobacterium
## 155 Oceanospirillaceae         Marinobacterium
## 163 Oceanospirillaceae         Marinobacterium
## 165 Oceanospirillaceae         Marinobacterium
## 166 Oceanospirillaceae         Marinobacterium
## 179     Cryomorphaceae                Wandonia
## 183     Cryomorphaceae                Wandonia
## 184     Cryomorphaceae                Wandonia
## 188     Cryomorphaceae                Wandonia
## 189     Cryomorphaceae                Wandonia
## 195     Cryomorphaceae                Wandonia
## 201     Cryomorphaceae                Wandonia
## 246  Flavobacteriaceae              Aquibacter
## 247  Flavobacteriaceae              Aquibacter
## 251  Flavobacteriaceae              Aquibacter
## 253  Flavobacteriaceae              Aquibacter
## 257  Flavobacteriaceae              Aquibacter
## 262  Flavobacteriaceae              Aquibacter
## 268  Flavobacteriaceae              Aquibacter
## 269  Flavobacteriaceae              Aquibacter
## 270  Flavobacteriaceae              Aquibacter
## 291 Oceanospirillaceae             Marinomonas
## 300 Oceanospirillaceae             Marinomonas
## 301 Oceanospirillaceae             Marinomonas
## 321   Alteromonadaceae         Aestuariibacter
## 333   Alteromonadaceae         Aestuariibacter
## 340   Alteromonadaceae         Aestuariibacter
## 376       Vibrionaceae                  Vibrio
## 388       Vibrionaceae                  Vibrio
## 395       Vibrionaceae                  Vibrio
## 396       Vibrionaceae                  Vibrio
## 397       Vibrionaceae                  Vibrio
## 399       Vibrionaceae                  Vibrio
## 400       Vibrionaceae                  Vibrio
## 401       Vibrionaceae                  Vibrio
## 402       Vibrionaceae                  Vibrio
## 403       Vibrionaceae                  Vibrio
## 405       Vibrionaceae                  Vibrio
## 407       Vibrionaceae                  Vibrio
## 408       Vibrionaceae                  Vibrio
```

``` r
cg_filter1=subset_samples(cg_filt,innoc=="innoculum" | strain_26 =="A1" | strain_26==
                            "A2" | strain_26== "B1" | strain_26 =="C1")
cg_filter1=subset_samples(cg_filt,strain_26 =="A1" | strain_26==
                            "A2" | strain_26== "B1" | strain_26 =="C1")

physeq_merged <- merge_samples(cg_filter1, group = "definition_26")
```

```
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
```

``` r
glom=tax_glom(physeq_merged, "species")

abund=transform_sample_counts(glom, function(x) (x/sum(x)))
keep_taxa <- taxa_sums(abund) / nsamples(abund) > 0.01
ps_filtered <- prune_taxa(keep_taxa, abund)
ps_filtered=ps_prune(abund, min.samples = 0, min.reads = 0, min.abundance = 0.01)
```

```
## 411 features grouped as 'Others' in the output
```

``` r
tax_table(ps_filtered)
```

```
## Taxonomy Table:     [15 taxa by 7 taxonomic ranks]:
##        Kingdom    phyla            class                 family             
## Otu6   "Bacteria" "Proteobacteria" "Alphaproteobacteria" "SAR11"            
## Otu19  "Bacteria" "Bacteroidetes"  "Flavobacteriia"      "Flavobacteriales" 
## Otu24  "Bacteria" "Bacteroidetes"  "Flavobacteriia"      "Flavobacteriales" 
## Otu2   "Bacteria" "Bacteroidetes"  "Flavobacteriia"      "Flavobacteriales" 
## Otu877 "Bacteria" "Bacteroidetes"  "Flavobacteriia"      "Flavobacteriales" 
## Otu10  "Bacteria" "Bacteroidetes"  "Flavobacteriia"      "Flavobacteriales" 
## Otu18  "Bacteria" "Bacteroidetes"  "Flavobacteriia"      "Flavobacteriales" 
## Otu11  "Bacteria" "Bacteroidetes"  "Flavobacteriia"      "Flavobacteriales" 
## Otu5   "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Alteromonadales"  
## Otu7   "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Vibrionales"      
## Otu4   "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Oceanospirillales"
## Otu12  "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Oceanospirillales"
## Otu14  "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Oceanospirillales"
## Otu16  "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Oceanospirillales"
## Others NA         NA               NA                    NA                 
##        genus                              species                   strain
## Otu6   "SAR11"                            "Candidatus Pelagibacter" NA    
## Otu19  "Cryomorphaceae"                   "Wandonia"                NA    
## Otu24  "Flavobacteriaceae"                "Aquibacter"              NA    
## Otu2   "Flavobacteriaceae"                "Tenacibaculum"           NA    
## Otu877 "Flavobacteriaceae"                "Polaribacter"            NA    
## Otu10  "Flavobacteriaceae"                "Cellulophaga"            NA    
## Otu18  "Cryomorphaceae"                   "Phaeocystidibacter"      NA    
## Otu11  "Cryomorphaceae"                   "Owenweeksia"             NA    
## Otu5   "Alteromonadaceae"                 "Aestuariibacter"         NA    
## Otu7   "Vibrionaceae"                     "Vibrio"                  NA    
## Otu4   "Oceanospirillaceae"               "Marinomonas"             NA    
## Otu12  "Oceanospirillaceae"               "Neptuniibacter"          NA    
## Otu14  "Oceanospirillaceae"               "Marinobacterium"         NA    
## Otu16  "Oceanospirillales_incertae_sedis" "Pseudohongiella"         NA    
## Others NA                                 NA                        NA
```

``` r
all_melt <- psmelt(ps_filtered)
summary(as.factor(all_melt$species))
```

```
##         Aestuariibacter              Aquibacter Candidatus Pelagibacter 
##                       2                       2                       2 
##            Cellulophaga         Marinobacterium             Marinomonas 
##                       2                       2                       2 
##          Neptuniibacter             Owenweeksia      Phaeocystidibacter 
##                       2                       2                       2 
##            Polaribacter         Pseudohongiella           Tenacibaculum 
##                       2                       2                       2 
##                  Vibrio                Wandonia                    NA's 
##                       2                       2                       2
```

``` r
length(summary(as.factor(all_melt$species)))
```

```
## [1] 15
```

``` r
all_melt <- all_melt %>%
  group_by(Sample, species) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")


summary(as.factor(all_melt$species))
```

```
##         Aestuariibacter              Aquibacter Candidatus Pelagibacter 
##                       2                       2                       2 
##            Cellulophaga         Marinobacterium             Marinomonas 
##                       2                       2                       2 
##          Neptuniibacter             Owenweeksia      Phaeocystidibacter 
##                       2                       2                       2 
##            Polaribacter         Pseudohongiella           Tenacibaculum 
##                       2                       2                       2 
##                  Vibrio                Wandonia                    NA's 
##                       2                       2                       2
```

``` r
length(summary(as.factor(all_melt$species)))
```

```
## [1] 15
```

``` r
all_melt <- all_melt %>%
  mutate(Sample = fct_relevel(Sample,
                              "Coculture", "Algal-cell associated"))

tol_muted_15 <- c(
  "#88CCEE", "#44AA99", "#117733", "#332288",
  "#DDCC77", "#999933", "#CC6677", "#882255",
  "#AA4499", "#6699CC", "#E69F00", "#D55E00",
  "#F0E442", "#CC79A7", "#000000"
)
tol_muted_15 <- c(
  "#88CCEE", "#44AA99", "#117733", "#332288",
  "#DDCC77", "#999933", "#CC6677", "#882255",
  "#AA4499", "#DDDDDD", "#6699CC", "#888888",
  "#E69F00", "#F0E442", "black"
)

tol_muted_14 <- c(
  "#88CCEE", "#44AA99", "#117733", "#332288",
  "#DDCC77", "#999933", "#CC6677", "#882255",
  "#AA4499", "#DDDDDD", "#6699CC", "#888888",
  "#E69F00", "black"
)

tol_muted_15 <- c(
  "#88CCEE", "#44AA99", "#117733", "#332288",
  "#DDCC77", "#999933", "#CC6677", "#882255",
  "#AA4499", "#DDDDDD", "#6699CC",
  "#E69F00", "#F0E442", "black"
)


all_melt$species <- as.factor(all_melt$species) %>%
  fct_relevel("Other", after = Inf)
```

```
## Warning: 1 unknown level in `f`: Other
```

``` r
all_figure=ggplot(all_melt, aes(x = Sample, y = Abundance, alluvium = species, stratum = species, fill = species)) +
  geom_flow(alpha = 0.7,na.rm=T) +
  geom_stratum() +
  theme_minimal(base_size=12) +
  theme(legend.position = 'bottom') +
  xlab("A1 sample flow") +
    scale_fill_manual(values = tol_muted_15) +
  ylab("Relative Abundance")
all_figure
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

``` r
summary(as.factor(all_melt$species))
```

```
##         Aestuariibacter              Aquibacter Candidatus Pelagibacter 
##                       2                       2                       2 
##            Cellulophaga         Marinobacterium             Marinomonas 
##                       2                       2                       2 
##          Neptuniibacter             Owenweeksia      Phaeocystidibacter 
##                       2                       2                       2 
##            Polaribacter         Pseudohongiella           Tenacibaculum 
##                       2                       2                       2 
##                  Vibrio                Wandonia                    NA's 
##                       2                       2                       2
```

``` r
length(summary(as.factor(all_melt$species)))
```

```
## [1] 15
```

``` r
library(ggalluvial)
library(ggplot2)


#A1 
cg_filter1=subset_samples(cg_filt,innoc=="innoculum" | strain_26 =="A1")
physeq_merged <- merge_samples(cg_filter1, group = "definition_26")
```

```
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
```

``` r
glom=tax_glom(physeq_merged, "species")
abund=transform_sample_counts(glom, function(x) (x/sum(x)))
a1_melt <- psmelt(abund)

keep_classes <- c("Aestuariibacter", "Aquibacter","Candidatus Pelagibacter", 
                  "Cellulophaga",
                  "Marinobacterium", "Marinomonas","Neptuniibacter", "Other",
                  "Owenweeksia", "Polaribacter","Pseudohongiella",
                  "Tenacibaculum", "Vibrio",
                  "Wandonia")
a1_summed <- a1_melt %>%
  mutate(species = fct_other(species, keep = keep_classes, other_level = "Other"))

a1_summed <- a1_summed %>%
  group_by(Sample, species) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")
library(forcats)
a1_summed <- a1_summed %>%
  mutate(Sample = fct_relevel(Sample,
                              "Innoculum", "Coculture", "Algal-cell associated"))
tol_muted <- c(
  "#88CCEE", "#44AA99", "#117733", "#332288",
  "#DDCC77", "#999933", "#CC6677", "#882255",
  "#AA4499", "#DDDDDD", "#6699CC", "#888888", 'black'
)

a1_figure=ggplot(a1_summed, aes(x = Sample, y = Abundance, alluvium = species, stratum = species, fill = species)) +
  geom_flow(alpha = 0.7,na.rm=T) +
  geom_stratum() +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'none') +
  xlab("A1") +
  scale_fill_manual(values = tol_muted_15) +
  ylab("Relative Abundance")
a1_figure
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-17-2.png)<!-- -->

``` r
a1_figure2=ggplot(a1_summed, aes(x = Sample, y = Abundance, alluvium = species, stratum = species, fill = species)) +
  geom_flow(alpha = 0.7,na.rm=T) +
  geom_stratum() +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  xlab("A1") +
  scale_fill_manual(values = tol_muted_15) +
  ylab("Relative Abundance")
a1_figure2
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-17-3.png)<!-- -->

``` r
#ggsave('flowplot_legend.svg', plot=a1_figure2)

#A2
cg_filter1=subset_samples(cg_filt,innoc=="innoculum" | strain_26 =="A2")
physeq_merged <- merge_samples(cg_filter1, group = "definition_26")
```

```
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
```

``` r
glom=tax_glom(physeq_merged, "species")
abund=transform_sample_counts(glom, function(x) (x/sum(x)))
a2_melt <- psmelt(abund)

a2_summed <- a2_melt %>%
  mutate(species = fct_other(species, keep = keep_classes, other_level = "Other"))

a2_summed <- a2_summed %>%
  group_by(Sample, species) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

library(forcats)
a2_summed <- a2_summed %>%
  mutate(Sample = fct_relevel(Sample,
                                     "Innoculum", "Coculture", "Algal-cell associated"))

a2_figure=ggplot(a2_summed, aes(x = Sample, y = Abundance, alluvium = species, stratum = species, fill = species)) +
  geom_flow(alpha = 0.7) +
  geom_stratum() +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'none') +
  xlab("A2") +
  scale_fill_manual(values = tol_muted_15) +
  ylab("Relative Abundance")
a2_figure
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-17-4.png)<!-- -->

``` r
#B1
cg_filter1=subset_samples(cg_filt,innoc=="innoculum" | strain_26 =="B1")
physeq_merged <- merge_samples(cg_filter1, group = "definition_26")
```

```
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
```

``` r
glom=tax_glom(physeq_merged, "species")
abund=transform_sample_counts(glom, function(x) (x/sum(x)))
b1_melt <- psmelt(abund)
b1_summed <- b1_melt %>%
  mutate(species = fct_other(species, keep = keep_classes, other_level = "Other"))

b1_summed <- b1_summed %>%
  group_by(Sample, species) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

b1_summed <- b1_summed %>%
  mutate(Sample = fct_relevel(Sample,
                              "Innoculum", "Coculture", "Algal-cell associated"))

b1_figure=ggplot(b1_summed, aes(x = Sample, y = Abundance, alluvium = species, stratum = species, fill = species)) +
  geom_flow(alpha = 0.7) +
  geom_stratum() +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'none') +
  xlab("B1") +
  scale_fill_manual(values = tol_muted_15) +
  ylab("Relative Abundance")
b1_figure
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-17-5.png)<!-- -->

``` r
# C1
cg_filter1=subset_samples(cg_filt,innoc=="innoculum" | strain_26 =="C1")
physeq_merged <- merge_samples(cg_filter1, group = "definition_26")
```

```
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
## Warning in asMethod(object): NAs introduced by coercion
```

``` r
glom=tax_glom(physeq_merged, "species")
abund=transform_sample_counts(glom, function(x) (x/sum(x)))

c1_melt <- psmelt(abund)
c1_summed <- c1_melt %>%
  mutate(species = fct_other(species, keep = keep_classes, other_level = "Other"))
c1_summed <- c1_summed %>%
  group_by(Sample, species) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

c1_summed <- c1_summed %>%
  mutate(Sample = fct_relevel(Sample,
                              "Innoculum", "Coculture", "Algal-cell associated"))

c1_figure=ggplot(c1_summed, aes(x = Sample, y = Abundance, alluvium = species, stratum = species, fill = species)) +
  geom_flow(alpha = 0.7) +
  geom_stratum() +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'none') +
  xlab("C1") +
  scale_fill_manual(values = tol_muted_15) +
  ylab("Relative Abundance")
c1_figure
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-17-6.png)<!-- -->

``` r
cowplot::plot_grid(a1_figure, a2_figure, b1_figure,c1_figure)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/unnamed-chunk-17-7.png)<!-- -->

``` r
cowww=cowplot::plot_grid(a1_figure, a2_figure, b1_figure,c1_figure,ncol =4)

#ggsave('ribbon_plot_genus.svg', plot=cowww, dpi=300, height=4, width=16)
#ggsave('ribbon_plot_genus_leg.svg', plot=a1_figure2, dpi=300)

#cowplot::plot_grid(a1_figure, a2_figure, b1_figure,c1_figure, 
  #                 ncol=1)
```
