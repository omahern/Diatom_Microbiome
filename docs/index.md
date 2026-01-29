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

```r
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
```


# Physiology
## Growth Rates
### Bacteria 

```r
read<-read.csv(file="datafiles/Bact_Growth_rates.csv",
               header=T,row.names=1)

boxplot(read$Growth~read$Geno)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/bact-1-1.png)<!-- -->

```r
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

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/bact-1-2.png)<!-- -->

```r
boxplot(read$Growth~read$Pop)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/bact-1-3.png)<!-- -->

```r
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

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/bact-1-4.png)<!-- -->

### Diatoms 

#### edited and used as of 8 jan

```r
read<-read.csv(file="datafiles/TR_Growth_rates_December25.csv",
               header=T,row.names=1)

subset=subset(read, treatment!="RD4_YE5" & treatment!="RA5_RD4" & treatment!="RD4_NB7")
subset=subset(read, treatment!="NB7_22" & treatment!="NB7_AX-22" & treatment!="RD4_22" & treatment!="RD4_AX-22")



growth_order=c("RA4_AX","RA4","RA5_AX","RA5","RD4_AX","RD4",
               'YC5_AX','YC5',"YE5_AX","YE5","NB4_AX","NB4","NB6_AX","NB6","NB7_AX",
               "NB7","RA5_RD4", "RD4_YE5", "RD4_NB7")

 subset$names2=paste(subset$isolate, subset$temp)

 
 
 
ggplot1=ggplot(subset, aes(x=isolate, y=growth_rate, 
                   fill=culture, 
                   color=population)) +
  geom_boxplot(alpha=.4) +
  theme_classic() + 
  scale_fill_manual(values=c('white','gray'))+
  scale_color_manual(values=c('black','navy','firebrick', 'forestgreen')) +
# scale_fill_manual(values=c("gray70",'navy',"gray70",'navy',"gray70",'navy',"gray70",'navy','gray70','firebrick','gray70','firebrick','gray70','forestgreen','gray70','forestgreen','gray70','forestgreen','gray70','forestgreen')) +
 # scale_color_manual(values=c('black','navy','black','navy','black','navy','black','navy','black','firebrick','black','firebrick','black','forestgreen','black','forestgreen','black','forestgreen','black','forestgreen' )) + 
  theme(plot.margin=unit(c(0,0,0,0.25),"cm"))+
   facet_grid(clonal~.,scales='free',space='free') +
  theme(legend.position="bottom",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 
 # scale_y_continuous(name='Growth Rate')+
  ylab(expression("Growth rate " * mu *  ~day^-1)) +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  #stat_compare_means(growth_rate~treatment, group.by='group',data=subset, 'anova')+
  scale_x_discrete(name='') + coord_flip()
ggplot1


compare_means(growth_rate~treatment, group.by='isolate2',data=subset, 'anova',paired=TRUE,p.adjust.method = 'none')
test=compare_means(growth_rate~treatment, group.by='isolate2',data=subset, 'anova',paired=TRUE,p.adjust.method = 'none')
test=compare_means(growth_rate~isolate, group.by='names',data=subset, 'anova',paired=TRUE,p.adjust.method = 'none')
test=compare_means(growth_rate~clonal, group.by='isolate2',data=subset, 'anova',paired=TRUE,p.adjust.method = 'none')
test=compare_means(growth_rate~clonal, data=subset, 'anova',paired=TRUE,p.adjust.method = 'none', group.by='isolate2')
test

ggplot1+stat_pvalue_manual(data=test,label='p.adj' )
```

### used 10 january


```r
read<-read.csv(file="datafiles/TR_Growth_rates_December25.csv",
               header=T,row.names=1)

subset=subset(read, treatment!="RD4_YE5" & treatment!="RA5_RD4" & treatment!="RD4_NB7")
subset=subset(read, treatment!="NB7_22" & treatment!="NB7_AX-22" & treatment!="RD4_22" & treatment!="RD4_AX-22")
subset=subset(read, isolate!="NB4" & isolate !="NB6" & isolate !="RA4" & isolate !="YC5")


growth_order=c("RA4_AX","RA4","RA5_AX","RA5","RD4_AX","RD4",
               'YC5_AX','YC5',"YE5_AX","YE5","NB4_AX","NB4","NB6_AX","NB6","NB7_AX",
               "NB7","RA5_RD4", "RD4_YE5", "RD4_NB7")

 subset$names2=paste(subset$isolate, subset$temp)
subset2=subset
 mono=subset(subset2, clonal=='monoclonal')
 
 
ggplot1=ggplot(subset, aes(x=isolate, y=growth_rate, 
                   fill=culture, 
                   color=population)) +
  geom_boxplot(alpha=.4) +
  theme_classic() + 
  scale_fill_manual(values=c('white','gray'))+
  scale_color_manual(values=c('black','navy','firebrick', 'forestgreen')) +
# scale_fill_manual(values=c("gray70",'navy',"gray70",'navy',"gray70",'navy',"gray70",'navy','gray70','firebrick','gray70','firebrick','gray70','forestgreen','gray70','forestgreen','gray70','forestgreen','gray70','forestgreen')) +
 # scale_color_manual(values=c('black','navy','black','navy','black','navy','black','navy','black','firebrick','black','firebrick','black','forestgreen','black','forestgreen','black','forestgreen','black','forestgreen' )) + 
  theme(plot.margin=unit(c(0,0,0,0.25),"cm"))+
   facet_grid(clonal~.,scales='free',space='free') +
 #theme(legend.position="bottom",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 
      theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 

 # scale_y_continuous(name='Growth Rate')+
  ylab(expression("Growth rate " * mu *  ~day^-1)) +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  #stat_compare_means(growth_rate~treatment, group.by='group',data=subset, 'anova')+
  scale_x_discrete(name='') +coord_flip()
ggplot1
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fcm-growth-2-1.png)<!-- -->

#### with fv.fm

```r
read=read.csv(file='/Users/oliviaahern/Documents/GitHub/MS2/ms2_fvfm.csv',header=T,row.names=1)
dim(read)
```

```
## [1] 48  5
```

```r
head(read)
```

```
##                X.1     fv strain  bact population
## RA4-A RA4-A_T5.000 0.5173    RA4 xenic       PopA
## RA4-B RA4-B_T5.000 0.5206    RA4 xenic       PopA
## RA4-C RA4-C_T5.000 0.5189    RA4 xenic       PopA
## RA5-A RA5-A_T5.000 0.5550    RA5 xenic       PopA
## RA5-B RA5-B_T5.000 0.5523    RA5 xenic       PopA
## RA5-C RA5-C_T5.000 0.5631    RA5 xenic       PopA
```

```r
read2=subset(read, strain!="NB4" & strain !="NB6" & strain !="RA4" & strain !="YC5")

gplot2=ggplot(read2, 
             aes(x=strain, y=fv,
                 color=factor(population), fill=factor(bact)))+
  geom_boxplot(alpha=.4) +
  theme_classic() + 
  scale_fill_manual(values=c('white','gray'))+
  scale_color_manual(values=c("navy","firebrick","forestgreen")) +
# theme(legend.position="bottom",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 
    theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 

 # scale_y_continuous(name='Growth Rate')+
  ylab("(Fm/Fo)/Fm") +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  #stat_compare_means(growth_rate~treatment, group.by='group',data=subset, 'anova')+
  scale_x_discrete(name='')+ coord_flip()
 none
```

```
## function (.x, .p, ...) 
## {
##     every(.x, negate(.p), ...)
## }
## <bytecode: 0x119b0e2e8>
## <environment: namespace:purrr>
```

```r
cowplot::plot_grid(ggplot1, gplot2, 
                   rel_heights = c(1.7,1),
                   nrow=2,
                   align='h',
                   axis='l',
                   byrow=TRUE, labels='auto')
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fcm-growth-a-1.png)<!-- -->

```r
       #           rel_widths = c(1.7,1))
```

##### stats of fv/fm

```r
a=aov(read$fv ~read$strain*read$bact)
tuk=TukeyHSD(a)
par(mar=c(5,10,1,1))
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fcm-growth-b-1.png)<!-- -->![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fcm-growth-b-2.png)<!-- -->![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fcm-growth-b-3.png)<!-- -->

```r
data=data.frame(read2)
data %>%
  group_by(strain, bact) %>%
  get_summary_stats(fv, type = "mean_sd")
```

```
## # A tibble: 8 × 6
##   strain bact   variable     n  mean    sd
##   <chr>  <chr>  <fct>    <dbl> <dbl> <dbl>
## 1 NB7    axenic fv           3 0.48  0.014
## 2 NB7    xenic  fv           3 0.502 0.008
## 3 RA5    axenic fv           3 0.545 0.003
## 4 RA5    xenic  fv           3 0.557 0.006
## 5 RD4    axenic fv           3 0.472 0.023
## 6 RD4    xenic  fv           3 0.425 0.009
## 7 YE5    axenic fv           3 0.475 0.006
## 8 YE5    xenic  fv           3 0.427 0.006
```

```r
data %>%
  group_by(strain, bact) %>%
  shapiro_test(fv)
```

```
## # A tibble: 8 × 5
##   strain bact   variable statistic     p
##   <chr>  <chr>  <chr>        <dbl> <dbl>
## 1 NB7    axenic fv           0.968 0.658
## 2 NB7    xenic  fv           0.979 0.725
## 3 RA5    axenic fv           0.967 0.650
## 4 RA5    xenic  fv           0.923 0.463
## 5 RD4    axenic fv           0.868 0.289
## 6 RD4    xenic  fv           0.996 0.887
## 7 YE5    axenic fv           0.941 0.530
## 8 YE5    xenic  fv           0.946 0.553
```

```r
model <- lm(fv ~ strain * bact, data = data)

data %>%
  group_by(strain) %>%
  anova_test(fv ~ bact, error = model)
```

```
## # A tibble: 4 × 8
##   strain Effect   DFn   DFd     F         p `p<.05`   ges
## * <chr>  <chr>  <dbl> <dbl> <dbl>     <dbl> <chr>   <dbl>
## 1 NB7    bact       1    16  5.65 0.03      "*"     0.261
## 2 RA5    bact       1    16  1.67 0.215     ""      0.094
## 3 RD4    bact       1    16 25.9  0.000109  "*"     0.618
## 4 YE5    bact       1    16 28.0  0.0000733 "*"     0.636
```

```r
mono=subset(subset2, clonal=='monoclonal' & temp =="18")
```


```r
a=aov(mono$growth~mono$population)
summary(a)
```

```
##                 Df Sum Sq Mean Sq F value   Pr(>F)    
## mono$population  2 0.8967  0.4483   45.91 2.15e-08 ***
## Residuals       21 0.2051  0.0098                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/stats-a-1.png)<!-- -->

```r
range(mono$growth_rate
      )
```

```
## [1] 0.7336728 1.4978085
```


```r
a=aov(subset$growth~subset$clonal)
summary(a)
```

```
##               Df Sum Sq Mean Sq F value Pr(>F)
## subset$clonal  1  0.127 0.12684   1.696    0.2
## Residuals     43  3.216 0.07478
```

```r
tuk=TukeyHSD(a)
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/stats-b-1.png)<!-- -->


```r
a=aov(subset$growth~subset$isolate*subset$culture)
summary(a)
tuk=TukeyHSD(a)
plot(tuk,las=2)

model <- lm(growth_rate ~ isolate * bact, data = subset2)

data %>%
  group_by(strain) %>%
  anova_test(fv ~ bact, error = model)
```


```r
par(mfrow=c(1,2))
boxplot(read$Growth~read$Geno)
a=aov(read$Growth~read$Geno)
summary(a)
tuk=TukeyHSD(a)
plot(tuk,las=2)
```


```r
boxplot(read$Growth~read$Pop)
a=aov(read$Growth~read$Pop)
summary(a)
tuk=TukeyHSD(a)
plot(tuk,las=2)
```







## fv/fm

```r
read=read.csv(file='datafiles/ms2_fvfm.csv',header=T,row.names=1)
dim(read)
```

```
## [1] 48  5
```

```r
head(read)
```

```
##                X.1     fv strain  bact population
## RA4-A RA4-A_T5.000 0.5173    RA4 xenic       PopA
## RA4-B RA4-B_T5.000 0.5206    RA4 xenic       PopA
## RA4-C RA4-C_T5.000 0.5189    RA4 xenic       PopA
## RA5-A RA5-A_T5.000 0.5550    RA5 xenic       PopA
## RA5-B RA5-B_T5.000 0.5523    RA5 xenic       PopA
## RA5-C RA5-C_T5.000 0.5631    RA5 xenic       PopA
```

```r
ggplot(read, 
             aes(x=strain, y=fv,
                 color=factor(population), fill=factor(bact)))+
  geom_boxplot(alpha=.4) +
  theme_classic() + 
  scale_fill_manual(values=c('white','gray'))+
  scale_color_manual(values=c("navy","firebrick","forestgreen")) +
  theme(legend.position="bottom",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size=10)) + 
 # scale_y_continuous(name='Growth Rate')+
  ylab("(Fm/Fo)/Fm") +
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))+
  #stat_compare_means(growth_rate~treatment, group.by='group',data=subset, 'anova')+
  scale_x_discrete(name='')  + coord_flip()
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-1-1.png)<!-- -->


```r
a=aov(read$fv ~read$strain*read$bact)
tuk=TukeyHSD(a)
par(mar=c(5,10,1,1))
plot(tuk,las=2)
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-2-1.png)<!-- -->![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-2-2.png)<!-- -->![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/fv-2-3.png)<!-- -->

```r
data=data.frame(read)
data %>%
  group_by(strain, bact) %>%
  get_summary_stats(fv, type = "mean_sd")
```

```
## # A tibble: 16 × 6
##    strain bact   variable     n  mean    sd
##    <chr>  <chr>  <fct>    <dbl> <dbl> <dbl>
##  1 NB4    axenic fv           3 0.525 0.018
##  2 NB4    xenic  fv           3 0.498 0.015
##  3 NB6    axenic fv           3 0.535 0.007
##  4 NB6    xenic  fv           3 0.51  0.006
##  5 NB7    axenic fv           3 0.48  0.014
##  6 NB7    xenic  fv           3 0.502 0.008
##  7 RA4    axenic fv           3 0.53  0.006
##  8 RA4    xenic  fv           3 0.519 0.002
##  9 RA5    axenic fv           3 0.545 0.003
## 10 RA5    xenic  fv           3 0.557 0.006
## 11 RD4    axenic fv           3 0.472 0.023
## 12 RD4    xenic  fv           3 0.425 0.009
## 13 YC5    axenic fv           3 0.563 0.008
## 14 YC5    xenic  fv           3 0.562 0.015
## 15 YE5    axenic fv           3 0.475 0.006
## 16 YE5    xenic  fv           3 0.427 0.006
```

```r
data %>%
  group_by(strain, bact) %>%
  shapiro_test(fv)
```

```
## # A tibble: 16 × 5
##    strain bact   variable statistic      p
##    <chr>  <chr>  <chr>        <dbl>  <dbl>
##  1 NB4    axenic fv           0.954 0.586 
##  2 NB4    xenic  fv           0.989 0.795 
##  3 NB6    axenic fv           0.884 0.337 
##  4 NB6    xenic  fv           0.979 0.720 
##  5 NB7    axenic fv           0.968 0.658 
##  6 NB7    xenic  fv           0.979 0.725 
##  7 RA4    axenic fv           0.822 0.168 
##  8 RA4    xenic  fv           1.00  0.967 
##  9 RA5    axenic fv           0.967 0.650 
## 10 RA5    xenic  fv           0.923 0.463 
## 11 RD4    axenic fv           0.868 0.289 
## 12 RD4    xenic  fv           0.996 0.887 
## 13 YC5    axenic fv           0.981 0.735 
## 14 YC5    xenic  fv           0.761 0.0255
## 15 YE5    axenic fv           0.941 0.530 
## 16 YE5    xenic  fv           0.946 0.553
```

```r
model <- lm(fv ~ strain * bact, data = data)

data %>%
  group_by(strain) %>%
  anova_test(fv ~ bact, error = model)
```

```
## # A tibble: 8 × 8
##   strain Effect   DFn   DFd      F          p `p<.05`      ges
## * <chr>  <chr>  <dbl> <dbl>  <dbl>      <dbl> <chr>      <dbl>
## 1 NB4    bact       1    32  9.58  0.004      "*"     0.23    
## 2 NB6    bact       1    32  7.46  0.01       "*"     0.189   
## 3 NB7    bact       1    32  5.80  0.022      "*"     0.154   
## 4 RA4    bact       1    32  1.39  0.247      ""      0.042   
## 5 RA5    bact       1    32  1.71  0.2        ""      0.051   
## 6 RD4    bact       1    32 26.6   0.0000125  "*"     0.454   
## 7 YC5    bact       1    32  0.016 0.901      ""      0.000496
## 8 YE5    bact       1    32 28.7   0.00000698 "*"     0.473
```





## Growth Curves

### Diatoms

```r
read=read.csv(file='datafiles/rotula_cellcounts.csv',header=T,row.names=1)
dim(read)
```

```
## [1] 69  7
```

```r
treatment=as.factor(read$treatment[1:48])

dim(read)
```

```
## [1] 69  7
```

```r
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

#### Plot all black and white

```r
{par(mfrow=c(3,5),mar=c(5,4,1,1))
  # NB4
  errbar(time,mean[,1],yplus=mean[,1]+sd[,1],yminus=mean[,1]-sd[,1],
         ylim=c(100,50000),log='y',type='b',
         ylab="Cell Abundance (cells/mL)",
         xlab="Time (days)",
         yaxt='n')
  errbar(time,mean[,2],yplus=mean[,2]+sd[,2],yminus=mean[,2]-sd[,2],col='gray70',
         type='b',add=T)
  text(1,200,"NB4 18",font=2)
  axis(2,at=c(100,1000,12500,50000))
  segments(0,100,0,60000,col='gray90',lty=2)
  segments(2,100,2,60000,col='gray90',lty=2)
  segments(3,100,3,60000,col='gray90',lty=2)
  segments(4,100,4,60000,col='gray90',lty=2)
  segments(5,100,5,60000,col='gray90',lty=2)
  
  
  #NB6
  errbar(time,mean[,3],yplus=mean[,3]+sd[,3],yminus=mean[,3]-sd[,3],
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  errbar(time,mean[,4],yplus=mean[,4]+sd[,4],yminus=mean[,4]-sd[,4],col='gray70',
         type='b',add=T)
  text(1,200,"NB6 18", font=2)
  axis(2,at=c(100,1000,12500,50000))
  segments(0,100,0,60000,col='gray90',lty=2)
  segments(2,100,2,60000,col='gray90',lty=2)
  segments(3,100,3,60000,col='gray90',lty=2)
  segments(4,100,4,60000,col='gray90',lty=2)
  segments(5,100,5,60000,col='gray90',lty=2)
  
  #NB7
  errbar(time,mean[,5],yplus=mean[,5]+sd[,5],yminus=mean[,5]-sd[,5],
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  errbar(time,mean[,6],yplus=mean[,6]+sd[,6],yminus=mean[,6]-sd[,6],col='gray70',
         type='b',add=T)
  text(1,200,"NB7 18", font=2)
  axis(2,at=c(100,1000,12500,50000))
  segments(0,100,0,60000,col='gray90',lty=2)
  segments(2,100,2,60000,col='gray90',lty=2)
  segments(3,100,3,60000,col='gray90',lty=2)
  segments(4,100,4,60000,col='gray90',lty=2)
  segments(5,100,5,60000,col='gray90',lty=2)
  
  #RA4
  errbar(time,mean[,7],yplus=mean[,7]+sd[,7],yminus=mean[,7]-sd[,7],
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  errbar(time,mean[,8],yplus=mean[,8]+sd[,8],yminus=mean[,8]-sd[,8],col='gray70',
         type='b',add=T)
  text(1,200,"RA4 18", font =2)
  axis(2,at=c(100,1000,12500,50000))
  segments(0,100,0,60000,col='gray90',lty=2)
  segments(2,100,2,60000,col='gray90',lty=2)
  segments(3,100,3,60000,col='gray90',lty=2)
  segments(4,100,4,60000,col='gray90',lty=2)
  segments(5,100,5,60000,col='gray90',lty=2)
  
  #RA5
  errbar(time,mean[,9],yplus=mean[,9]+sd[,9],yminus=mean[,9]-sd[,9],
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  errbar(time,mean[,10],yplus=mean[,10]+sd[,10],yminus=mean[,10]-sd[,10],col='gray70',
         type='b',add=T)
  text(1,200,"RA5 18",font=2)
  axis(2,at=c(100,1000,12500,50000))
  segments(0,100,0,60000,col='gray90',lty=2)
  segments(2,100,2,60000,col='gray90',lty=2)
  segments(3,100,3,60000,col='gray90',lty=2)
  segments(4,100,4,60000,col='gray90',lty=2)
  segments(5,100,5,60000,col='gray90',lty=2)
  
  #RD4
  errbar(time,mean[,11],yplus=mean[,11]+sd[,11],yminus=mean[,11]-sd[,11],
         ylim=c(100,50000),log='y',type='b',
         ylab="Cell Abundance (cells/mL)",
         xlab="Time (days)",
         yaxt='n')
  errbar(time,mean[,12],yplus=mean[,12]+sd[,12],yminus=mean[,12]-sd[,12],col='gray70',
         type='b',add=T)
  text(1,200,"RD4 18",font=2)
  axis(2,at=c(100,1000,12500,50000))
  segments(0,100,0,60000,col='gray90',lty=2)
  segments(2,100,2,60000,col='gray90',lty=2)
  segments(3,100,3,60000,col='gray90',lty=2)
  segments(4,100,4,60000,col='gray90',lty=2)
  segments(5,100,5,60000,col='gray90',lty=2)
  
  #YC5
  errbar(time,mean[,13],yplus=mean[,13]+sd[,13],yminus=mean[,13]-sd[,13],
         ylim=c(10,10000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  errbar(time,mean[,14],yplus=mean[,14]+sd[,14],yminus=mean[,14]-sd[,14],col='gray70',
         type='b',add=T)
  text(1,20,"YC5 18",font=2)
  axis(2,at=c(10,100,1250,5000))
  segments(0,1,0,60000,col='gray90',lty=2)
  segments(2,1,2,60000,col='gray90',lty=2)
  segments(3,1,3,60000,col='gray90',lty=2)
  segments(4,1,4,60000,col='gray90',lty=2)
  segments(5,1,5,60000,col='gray90',lty=2)
  
  #YE5
  errbar(time,mean[,15],yplus=mean[,15]+sd[,15],yminus=mean[,15]-sd[,15],
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  errbar(time,mean[,16],yplus=mean[,16]+sd[,16],yminus=mean[,16]-sd[,16],col='gray70',
         type='b',add=T)
  text(1,200,"YE5 18*", font=2)
  axis(2,at=c(100,1000,12500,50000))
  segments(0,100,0,60000,col='gray90',lty=2)
  segments(2,100,2,60000,col='gray90',lty=2)
  segments(3,100,3,60000,col='gray90',lty=2)
  segments(4,100,4,60000,col='gray90',lty=2)
  segments(5,100,5,60000,col='gray90',lty=2)
  
  

  
  # RA4-RD4
  errbar(time,mean_db[,1],yplus=mean_db[,1]+sd_db[,1],yminus=mean_db[,1]-sd_db[,1],
         ylim=c(100,60000),log='y',type='b',
         ylab="Cell Abundance (cells/mL)",
         xlab="Time (days)",
         yaxt='n')
  text(1,200,"RD4 + RA4",font=2)
  axis(2,at=c(100,1000,12500,50000))
  segments(0,100,0,60000,col='gray90',lty=2)
  segments(2,100,2,60000,col='gray90',lty=2)
  segments(3,100,3,60000,col='gray90',lty=2)
  segments(4,100,4,60000,col='gray90',lty=2)
  segments(5,100,5,60000,col='gray90',lty=2)
  
  # RD4-NB7
  errbar(time,mean_db[,2],yplus=mean_db[,2]+sd_db[,2],yminus=mean_db[,2]-sd_db[,2],
         ylim=c(100,60000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  text(1,200,"RD4 + NB7",font=2)
  axis(2,at=c(100,1000,12500,50000))
  segments(0,100,0,60000,col='gray90',lty=2)
  segments(2,100,2,60000,col='gray90',lty=2)
  segments(3,100,3,60000,col='gray90',lty=2)
  segments(4,100,4,60000,col='gray90',lty=2)
  segments(5,100,5,60000,col='gray90',lty=2)
  
  # YD4_YE5
  errbar(time,mean_db[,3],yplus=mean_db[,3]+sd_db[,3],yminus=mean_db[,3]-sd_db[,3],
         ylim=c(100,60000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  text(1,200,"RD4 + YE5",font=2)
  axis(2,at=c(100,1000,12500,50000))
  segments(0,100,0,60000,col='gray90',lty=2)
  segments(2,100,2,60000,col='gray90',lty=2)
  segments(3,100,3,60000,col='gray90',lty=2)
  segments(4,100,4,60000,col='gray90',lty=2)
  segments(5,100,5,60000,col='gray90',lty=2)
}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/curves-2-1.png)<!-- -->

#### Plot all colored

```r
read=read.csv(file='datafiles/rotula_cellcounts.csv',header=T,row.names=1)
dim(read)
```

```
## [1] 69  7
```

```r
treatment=as.factor(read$treatment[1:48])

dim(read)
```

```
## [1] 69  7
```

```r
growth_data=data.frame(read[1:48,3:7])

# get mean and standard deviation
mean=aggregate(growth_data,by=list(treatment),FUN=mean,na.rm=T)
mean=t(mean[,2:6])
sd=aggregate(growth_data,by=list(treatment),FUN=sd,na.rm=T)
sd=t(sd[,2:6])
# set time
time=c(0,2,3,4,5)

{par(mfrow=c(4,3),mar=c(2,2,1,1))
  #RA4
  errbar(time,mean[,7],yplus=mean[,7]+sd[,7],yminus=mean[,7]-sd[,7],
         ylim=c(100,50000),log='y',type='b', col='navy', pch=21, bg='gray',
         ylab="Cell Abundance (cells/mL)",
         errbar.col='navy', cex=1.3,xaxt='n',
         xlab=" ",
         yaxt='n')
  errbar(time,mean[,8],yplus=mean[,8]+sd[,8],yminus=mean[,8]-sd[,8],col='navy', pch=21, bg='white',
         type='b',add=T,errbar.col='navy', cex=1.3)
  text(1,200,"RA4", font =2)
  axis(2,at=c(100,1000,12500,50000))
  axis(1, at=c(0:5),label=FALSE)
  legend("bottomright", legend=c("xenic",'axenic'), pch=21, pt.bg=c('grey','white'), bty='n')
  
  #RA5
  errbar(time,mean[,9],yplus=mean[,9]+sd[,9],yminus=mean[,9]-sd[,9],
         ylim=c(100,50000),log='y',type='b',col='navy', pch=21, bg='gray',errbar.col = 'navy',cex=1.3,
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,10],yplus=mean[,10]+sd[,10],yminus=mean[,10]-sd[,10],col='navy', pch=21, bg='white',errbar.col = 'navy',cex=1.3,
         type='b',add=T)
  text(1,200,"RA5",font=2)
  axis(2,at=c(100,1000,12500,50000))
  axis(1, at=c(0:5),label=FALSE)

  
  #RD4
  errbar(time,mean[,11],yplus=mean[,11]+sd[,11],yminus=mean[,11]-sd[,11],
         ylim=c(100,50000),log='y',type='b',col='navy', pch=21, bg='gray',errbar.col = 'navy',cex=1.3,
         ylab="Cell Abundance (cells/mL)",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,12],yplus=mean[,12]+sd[,12],yminus=mean[,12]-sd[,12],col='navy', pch=21,cex=1.3, bg='white',errbar.col = 'navy',
         type='b',add=T)
  text(1,200,"RD4",font=2)
  axis(2,at=c(100,1000,12500,50000))
    axis(1, at=c(0:5),label=FALSE)

  
   

  
 

  
  #YC5
  errbar(time,mean[,13],yplus=mean[,13]+sd[,13],yminus=mean[,13]-sd[,13],col='firebrick', pch=21, bg='gray',errbar.col = 'firebrick',cex=1.3,
         ylim=c(10,10000),log='y',type='b',
         ylab="Cell Abundance (cells/mL)",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,14],yplus=mean[,14]+sd[,14],yminus=mean[,14]-sd[,14],col='firebrick', pch=21, bg='white',errbar.col = 'firebrick',cex=1.3,
         type='b',add=T)
  text(1,20,"YC5",font=2)
  axis(2,at=c(10,100,1250,5000))
    axis(1, at=c(0:5),label=FALSE)

  
  #YE5
  errbar(time,mean[,15],yplus=mean[,15]+sd[,15],yminus=mean[,15]-sd[,15],col='firebrick', pch=21, bg='gray',errbar.col = 'firebrick',cex=1.3,
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,16],yplus=mean[,16]+sd[,16],yminus=mean[,16]-sd[,16],col='firebrick', pch=21, bg='white',errbar.col = 'firebrick',cex=1.3,
         type='b',add=T)
  text(1,200,"YE5", font=2)
  axis(2,at=c(100,1000,12500,50000))
    axis(1, at=c(0:5),label=FALSE)
    
    plot.new()

  # NB4
  errbar(time,mean[,1],yplus=mean[,1]+sd[,1],yminus=mean[,1]-sd[,1],col='forestgreen', pch=21, bg='gray',errbar.col = 'forestgreen',cex=1.3,
         ylim=c(100,50000),log='y',type='b',
         ylab="Cell Abundance (cells/mL)",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,2],yplus=mean[,2]+sd[,2],yminus=mean[,2]-sd[,2],col='forestgreen', pch=21, bg='white',errbar.col = 'forestgreen',cex=1.3,
         type='b',add=T)
  text(1,200,"NB4",font=2)
  axis(2,at=c(100,1000,12500,50000))
    axis(1, at=c(0:5),label=FALSE)

  
  
  #NB6
  errbar(time,mean[,3],yplus=mean[,3]+sd[,3],yminus=mean[,3]-sd[,3],col='forestgreen', pch=21, bg='gray',errbar.col = 'forestgreen',cex=1.3,
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,4],yplus=mean[,4]+sd[,4],yminus=mean[,4]-sd[,4],col='forestgreen', pch=21, bg='white',errbar.col = 'forestgreen',cex=1.3,
         type='b',add=T)
  text(1,200,"NB6", font=2)
  axis(2,at=c(100,1000,12500,50000))
    axis(1, at=c(0:5),label=FALSE)

  
  #NB7
  errbar(time,mean[,5],yplus=mean[,5]+sd[,5],yminus=mean[,5]-sd[,5],col='forestgreen', pch=21, bg='gray',errbar.col = 'forestgreen',cex=1.3,
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,6],yplus=mean[,6]+sd[,6],yminus=mean[,6]-sd[,6],col='forestgreen', pch=21, bg='white',errbar.col = 'forestgreen',cex=1.3,
         type='b',add=T)
  text(1,200,"NB7", font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

  

  
  # RA4-RD4
  errbar(time,mean_db[,1],yplus=mean_db[,1]+sd_db[,1],yminus=mean_db[,1]-sd_db[,1],col='black', pch=21, bg='gray',errbar.col = 'black',cex=1.3,
         ylim=c(100,60000),log='y',type='b',
         ylab="Cell Abundance (cells/mL)",
         xlab=" ",
         yaxt='n')
  text(1,200,"RD4 + RA4",font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

  
  # RD4-NB7
  errbar(time,mean_db[,2],yplus=mean_db[,2]+sd_db[,2],yminus=mean_db[,2]-sd_db[,2],col='black', pch=21, bg='gray',errbar.col = 'black',cex=1.3,
         ylim=c(100,60000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  text(1,200,"RD4 + NB7",font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

  
  # YD4_YE5
  errbar(time,mean_db[,3],yplus=mean_db[,3]+sd_db[,3],yminus=mean_db[,3]-sd_db[,3],col='black', pch=21, bg='gray',errbar.col = 'black',cex=1.3,
         ylim=c(100,60000),log='y',type='b',
         ylab=" ",
         xlab=" ",
         yaxt='n')
  text(1,200,"RD4 + YE5",font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/curves-3-1.png)<!-- -->

#### Plot all colored 2

```r
read=read.csv(file='datafiles/rotula_cellcounts.csv',header=T,row.names=1)
dim(read)
```

```
## [1] 69  7
```

```r
treatment=as.factor(read$treatment[1:48])

dim(read)
```

```
## [1] 69  7
```

```r
growth_data=data.frame(read[1:48,3:7])

# get mean and standard deviation
mean=aggregate(growth_data,by=list(treatment),FUN=mean,na.rm=T)
mean=t(mean[,2:6])
sd=aggregate(growth_data,by=list(treatment),FUN=sd,na.rm=T)
sd=t(sd[,2:6])
# set time
time=c(0,2,3,4,5)

{par(mfrow=c(3,4),
     mar=c(2,2,1,1))
  #RA4
  errbar(time,mean[,7],yplus=mean[,7]+sd[,7],yminus=mean[,7]-sd[,7],
         ylim=c(100,50000),log='y',type='b', col='navy', pch=21, bg='gray',
         ylab="Cell Abundance (cells/mL)",
         errbar.col='navy', cex=1.3,xaxt='n',
         xlab=" ",
         yaxt='n')
  errbar(time,mean[,8],yplus=mean[,8]+sd[,8],yminus=mean[,8]-sd[,8],col='navy', pch=21, bg='white',
         type='b',add=T,errbar.col='navy', cex=1.3)
  text(1,200,"RA4", font =2)
  axis(2,at=c(100,1000,12500,50000))
  axis(1, at=c(0:5),label=FALSE)
  legend("bottomright", legend=c("xenic",'axenic'), pch=21, pt.bg=c('grey','white'), bty='n')
  
  #RA5
  errbar(time,mean[,9],yplus=mean[,9]+sd[,9],yminus=mean[,9]-sd[,9],
         ylim=c(100,50000),log='y',type='b',col='navy', pch=21, bg='gray',errbar.col = 'navy',cex=1.3,
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,10],yplus=mean[,10]+sd[,10],yminus=mean[,10]-sd[,10],col='navy', pch=21, bg='white',errbar.col = 'navy',cex=1.3,
         type='b',add=T)
  text(1,200,"RA5",font=2)
  axis(2,at=c(100,1000,12500,50000))
  axis(1, at=c(0:5),label=FALSE)

  
  #RD4
  errbar(time,mean[,11],yplus=mean[,11]+sd[,11],yminus=mean[,11]-sd[,11],
         ylim=c(100,50000),log='y',type='b',col='navy', pch=21, bg='gray',errbar.col = 'navy',cex=1.3,
         ylab="Cell Abundance (cells/mL)",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,12],yplus=mean[,12]+sd[,12],yminus=mean[,12]-sd[,12],col='navy', pch=21,cex=1.3, bg='white',errbar.col = 'navy',
         type='b',add=T)
  text(1,200,"RD4",font=2)
  axis(2,at=c(100,1000,12500,50000))
    axis(1, at=c(0:5),label=FALSE)

  
   

  
 

  
  #YC5
  errbar(time,mean[,13],yplus=mean[,13]+sd[,13],yminus=mean[,13]-sd[,13],col='firebrick', pch=21, bg='gray',errbar.col = 'firebrick',cex=1.3,
         ylim=c(10,10000),log='y',type='b',
         ylab="Cell Abundance (cells/mL)",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,14],yplus=mean[,14]+sd[,14],yminus=mean[,14]-sd[,14],col='firebrick', pch=21, bg='white',errbar.col = 'firebrick',cex=1.3,
         type='b',add=T)
  text(1,20,"YC5",font=2)
  axis(2,at=c(10,100,1250,5000))
    axis(1, at=c(0:5),label=FALSE)

  
  #YE5
  errbar(time,mean[,15],yplus=mean[,15]+sd[,15],yminus=mean[,15]-sd[,15],col='firebrick', pch=21, bg='gray',errbar.col = 'firebrick',cex=1.3,
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,16],yplus=mean[,16]+sd[,16],yminus=mean[,16]-sd[,16],col='firebrick', pch=21, bg='white',errbar.col = 'firebrick',cex=1.3,
         type='b',add=T)
  text(1,200,"YE5", font=2)
  axis(2,at=c(100,1000,12500,50000))
    axis(1, at=c(0:5),label=FALSE)
    
    plot.new()

  # NB4
  errbar(time,mean[,1],yplus=mean[,1]+sd[,1],yminus=mean[,1]-sd[,1],col='forestgreen', pch=21, bg='gray',errbar.col = 'forestgreen',cex=1.3,
         ylim=c(100,50000),log='y',type='b',
         ylab="Cell Abundance (cells/mL)",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,2],yplus=mean[,2]+sd[,2],yminus=mean[,2]-sd[,2],col='forestgreen', pch=21, bg='white',errbar.col = 'forestgreen',cex=1.3,
         type='b',add=T)
  text(1,200,"NB4",font=2)
  axis(2,at=c(100,1000,12500,50000))
    axis(1, at=c(0:5),label=FALSE)

  
  
  #NB6
  errbar(time,mean[,3],yplus=mean[,3]+sd[,3],yminus=mean[,3]-sd[,3],col='forestgreen', pch=21, bg='gray',errbar.col = 'forestgreen',cex=1.3,
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,4],yplus=mean[,4]+sd[,4],yminus=mean[,4]-sd[,4],col='forestgreen', pch=21, bg='white',errbar.col = 'forestgreen',cex=1.3,
         type='b',add=T)
  text(1,200,"NB6", font=2)
  axis(2,at=c(100,1000,12500,50000))
    axis(1, at=c(0:5),label=FALSE)

  
  #NB7
  errbar(time,mean[,5],yplus=mean[,5]+sd[,5],yminus=mean[,5]-sd[,5],col='forestgreen', pch=21, bg='gray',errbar.col = 'forestgreen',cex=1.3,
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,6],yplus=mean[,6]+sd[,6],yminus=mean[,6]-sd[,6],col='forestgreen', pch=21, bg='white',errbar.col = 'forestgreen',cex=1.3,
         type='b',add=T)
  text(1,200,"NB7", font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

  

  
  # RA4-RD4
  errbar(time,mean_db[,1],yplus=mean_db[,1]+sd_db[,1],yminus=mean_db[,1]-sd_db[,1],col='black', pch=21, bg='gray',errbar.col = 'black',cex=1.3,
         ylim=c(100,60000),log='y',type='b',
         ylab="Cell Abundance (cells/mL)",
         xlab=" ",
         yaxt='n')
  text(1,200,"RD4 + RA4",font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

  
  # RD4-NB7
  errbar(time,mean_db[,2],yplus=mean_db[,2]+sd_db[,2],yminus=mean_db[,2]-sd_db[,2],col='black', pch=21, bg='gray',errbar.col = 'black',cex=1.3,
         ylim=c(100,60000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  text(1,200,"RD4 + NB7",font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

  
  # YD4_YE5
  errbar(time,mean_db[,3],yplus=mean_db[,3]+sd_db[,3],yminus=mean_db[,3]-sd_db[,3],col='black', pch=21, bg='gray',errbar.col = 'black',cex=1.3,
         ylim=c(100,60000),log='y',type='b',
         ylab=" ",
         xlab=" ",
         yaxt='n')
  text(1,200,"RD4 + YE5",font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/curves-4-1.png)<!-- -->

#### Plot Subset of Strains (used)

```r
read=read.csv(file='datafiles/rotula_cellcounts.csv',header=T,row.names=1)
# dim(read)
treatment=as.factor(read$treatment[1:48])

# dim(read)
growth_data=data.frame(read[1:48,3:7])

# get mean and standard deviation
mean=aggregate(growth_data,by=list(treatment),FUN=mean,na.rm=T)
mean=t(mean[,2:6])
sd=aggregate(growth_data,by=list(treatment),FUN=sd,na.rm=T)
sd=t(sd[,2:6])
# set time
time=c(0,2,3,4,5)

{par(mfrow=c(2,4),
     mar=c(2,2,1,1))

  #RA5
  errbar(time,mean[,9],yplus=mean[,9]+sd[,9],yminus=mean[,9]-sd[,9],
         ylim=c(100,50000),log='y',type='b',col='navy', pch=21, bg='gray',errbar.col = 'navy',cex=1.3,
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,10],yplus=mean[,10]+sd[,10],yminus=mean[,10]-sd[,10],col='navy', pch=21, bg='white',errbar.col = 'navy',cex=1.3,
         type='b',add=T)
  text(1,200,"RA5",font=2)
  axis(2,at=c(100,1000,12500,50000))
  axis(1, at=c(0:5),label=FALSE)
  legend("bottomright", legend=c("xenic",'axenic'), pch=21, pt.bg=c('grey','white'), bty='n')

  
  #RD4
  errbar(time,mean[,11],yplus=mean[,11]+sd[,11],yminus=mean[,11]-sd[,11],
         ylim=c(100,50000),log='y',type='b',col='navy', pch=21, bg='gray',errbar.col = 'navy',cex=1.3,
         ylab="Cell Abundance (cells/mL)",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,12],yplus=mean[,12]+sd[,12],yminus=mean[,12]-sd[,12],col='navy', pch=21,cex=1.3, bg='white',errbar.col = 'navy',
         type='b',add=T)
  text(1,200,"RD4",font=2)
  axis(2,at=c(100,1000,12500,50000))
    axis(1, at=c(0:5),label=FALSE)



  
 

  
  #YE5
  errbar(time,mean[,15],yplus=mean[,15]+sd[,15],yminus=mean[,15]-sd[,15],col='firebrick', pch=21, bg='gray',errbar.col = 'firebrick',cex=1.3,
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,16],yplus=mean[,16]+sd[,16],yminus=mean[,16]-sd[,16],col='firebrick', pch=21, bg='white',errbar.col = 'firebrick',cex=1.3,
         type='b',add=T)
  text(1,200,"YE5", font=2)
  axis(2,at=c(100,1000,12500,50000))
    axis(1, at=c(0:5),label=FALSE)
    
 
  
  #NB7
  errbar(time,mean[,5],yplus=mean[,5]+sd[,5],yminus=mean[,5]-sd[,5],col='forestgreen', pch=21, bg='gray',errbar.col = 'forestgreen',cex=1.3,
         ylim=c(100,50000),log='y',type='b',
         ylab=" ",
         xlab=" ",xaxt='n',
         yaxt='n')
  errbar(time,mean[,6],yplus=mean[,6]+sd[,6],yminus=mean[,6]-sd[,6],col='forestgreen', pch=21, bg='white',errbar.col = 'forestgreen',cex=1.3,
         type='b',add=T)
  text(1,200,"NB7", font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)



  
  # RA4-RD4
  errbar(time,mean_db[,1],yplus=mean_db[,1]+sd_db[,1],yminus=mean_db[,1]-sd_db[,1],col='black', pch=21, bg='gray',errbar.col = 'black',cex=1.3,
         ylim=c(100,60000),log='y',type='b',
         ylab="Cell Abundance (cells/mL)",
         xlab=" ",
         yaxt='n')
  text(1,200,"RD4 + RA4",font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

  
  # RD4-NB7
  errbar(time,mean_db[,2],yplus=mean_db[,2]+sd_db[,2],yminus=mean_db[,2]-sd_db[,2],col='black', pch=21, bg='gray',errbar.col = 'black',cex=1.3,
         ylim=c(100,60000),log='y',type='b',
         ylab=" ",
         xlab="Time (days)",
         yaxt='n')
  text(1,200,"RD4 + NB7",font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

  
  # YD4_YE5
  errbar(time,mean_db[,3],yplus=mean_db[,3]+sd_db[,3],yminus=mean_db[,3]-sd_db[,3],col='black', pch=21, bg='gray',errbar.col = 'black',cex=1.3,
         ylim=c(100,60000),log='y',type='b',
         ylab=" ",
         xlab=" ",
         yaxt='n')
  text(1,200,"RD4 + YE5",font=2)
  axis(2,at=c(100,1000,12500,50000))
      axis(1, at=c(0:5),label=FALSE)

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/curves-5-1.png)<!-- -->

### bacteria (need to edit)
 
```
read=read.csv(file="/Users/ahern/Desktop/bact_growth.csv",header=T,row.names=1)

dim(read)

range(read$X5.00[8:43],na.rm=T)
diatoms=read$X5.00[8:43]
range(read$X0.00,na.rm=T)
mean(read$X0.00,na.rm=T)
sd(read$X0.00,na.rm=T)

diatoms=read$X5.00[8:46]
mean(diatoms,na.rm=T) # 976944.7
sd(diatoms,na.rm=T) # 1303581

ye5=read$X5.00[44:46]
mean(618605.0, 218731.7,530102.4)
s=c(618605.0, 218731.7,530102.4)
sd(s) # 618,605 +/- 210,033


bact=read$X5.00[1:7]
mean(bact)
sd(bact) # 382507.6 +/- 87,390.54

data=data.frame(read[,3:7])
plot(data[1,])
treatment=read$treat
mean1=aggregate(data,by=list(treatment),FUN=mean,na.rm=T)
mean=t(mean1[,2:6])
sd1=aggregate(data,by=list(treatment),FUN=sd,na.rm=T)
sd=t(sd1[,2:6])
# set time
time=c(0,2,3,4,5)

range(data,na.rm=T)

library(Hmisc)
sd=t(sd1[,2:6])

par(mfrow=c(2,2),mar=c(5,5,1,1))
  # Bact
  
  # NB
  errbar(time,mean[,1],yplus=mean[,1]+sd[,1],yminus=mean[,1]-sd[,1],
         type='b',ylim=c(1,4500), log='y',cex=1.5,
         ylab="Bacteria:Diatom ratio", pch=21, bg='forestgreen',
         xlab="Time (days)",cex.lab=1.2) 
  errbar(time,mean[,2],yplus=mean[,2]+sd[,2],yminus=mean[,2]-sd[,2],bg='forestgreen', pch =22,cex=1.5,
         type='b',add=T)
  errbar(time,mean[,3],yplus=mean[,3]+sd[,3],yminus=mean[,3]-sd[,3],bg='forestgreen', pch=23,cex=1.5,
         type='b',add=T)
  errbar(time,mean[,4],yplus=mean[,4]+sd[,4],yminus=mean[,4]-sd[,4],bg='green', pch=24,cex=1.5,
         type='b',add=T)
  legend('bottomleft',legend=c('NB4','NB6','NB7','NB7-22'), 
         pt.bg =c('forestgreen','forestgreen','forestgreen','green'),
         pch=c(21,22,23,24),bty='n')
  #text(4,200,"NB4")
  
  #
  # R
  errbar(time,mean[,5],yplus=mean[,5]+sd[,5],yminus=mean[,5]-sd[,5],
         log='y',type='b',ylim=c(1,4500), pch = 21, bg='navy',
         ylab="Bacteria:Diatom Ratio", cex.lab=1.2,cex=1.5,
         xlab="Time (days)")
  errbar(time,mean[,6],yplus=mean[,6]+sd[,6],yminus=mean[,6]-sd[,6],bg='navy', pch =22,cex=1.5,
         type='b',add=T)
  errbar(time,mean[,8],yplus=mean[,8]+sd[,8],yminus=mean[,8]-sd[,8],bg='navy', pch =23,cex=1.5,
         type='b',add=T)
  errbar(time,mean[,9],yplus=mean[,9]+sd[,9],yminus=mean[,9]-sd[,9],bg='blue', pch=24,cex=1.5,
         type='b',add=T)
  legend('bottomleft',legend=c('RA4','RA5','RD4','RD4-22'), 
         pt.bg =c('navy','navy','navy','blue'),
         pch=c(21,22,23,24),bty='n')
 # text(4,200,"NB4")
  
  # R
  errbar(time,mean[,12],yplus=mean[,12]+sd[,12],yminus=mean[,12]-sd[,12],
         log='y',type='b',ylim=c(1,50000), pch =21, bg='firebrick',
         ylab="Bacteria:Diatom Ratio", cex.lab=1.2,cex=1.5,
         xlab="Time (days)")
  errbar(time,mean[,13],yplus=mean[,13]+sd[,13],yminus=mean[,13]-sd[,13],bg='firebrick', pch = 22,cex=1.5,
         type='b',add=T)
  legend('bottomleft',legend=c('YC5','YE5'), 
         pt.bg =c('firebrick','firebrick'),
         pch=c(21,22),bty='n')
  
  # R
  errbar(time,mean[,7],yplus=mean[,7]+sd[,7],yminus=mean[,7]-sd[,7],
         log='y',type='b',ylim=c(1,4500), pch =24,cex=1.5,
         ylab="Diatom:Bacteria Ratio", bg= 'blue', cex.lab=1.2,
         xlab="Time (days)")
  errbar(time,mean[,10],yplus=mean[,10]+sd[,10],yminus=mean[,10]-sd[10],bg='green', pch = 24,cex=1.5,
         type='b',add=T)
  errbar(time,mean[,11],yplus=mean[,11]+sd[,11],yminus=mean[,11]-sd[11],bg='orange', pch = 24,cex=1.5,
         type='b',add=T)
  legend('bottomleft',legend=c('RD4 + RA5','RD4 + NB7', 'RD4 + YE5'), 
         pt.bg =c('blue','green','orange'),
         pch=24,bty='n')
  
dev.off()
  read$gen_time1
  

  
  
  ## try diatom bacteria ratio 
  ###
  read=read.csv(file='/Users/ahern/R/CHP_3/diatom_bacteria_ratio.csv',header=T,row.names=1)
  data=data.frame(read[,3:7])
  plot(data[1,])
  treatment=read$treat
  mean1=aggregate(data,by=list(treatment),FUN=mean,na.rm=T)
  mean=t(mean1[,2:6])
  sd1=aggregate(data,by=list(treatment),FUN=sd,na.rm=T)
  sd=t(sd1[,2:6])
  # set time
  time=c(0,2,3,4,5)
  
  time=c(0,2,3,4,5)
  
  library(Hmisc)
  
  par(mfrow=c(2,3),mar=c(5,5,1,1))
  # Bact
  errbar(time,mean[,1],yplus=mean[,1]+sd[,1],yminus=mean[,1]-sd[,1],
         log='y',type='b',ylim=c(119908.6,7787004.6),
         ylab="Cell Abundance (cells/mL)",
         xlab="Time (days)")
  errbar(time,mean[,2],yplus=mean[,2]+sd[,2],yminus=mean[,2]-sd[,2],col='gray70',
         type='b',add=T)
  text(1,5000000,"Bact 18C")
  text(1,4000000,"Bact 22C",col='gray70')
  
  # NB
  errbar(time,mean[,3],yplus=mean[,3]+sd[,3],yminus=mean[,3]-sd[,3],
         log='y',type='b',ylim=c(119908.6,7787004.6),
         ylab="Cell Abundance (cells/mL)", pch=21, bg='forestgreen',
         xlab="Time (days)") 
  errbar(time,mean[,4],yplus=mean[,4]+sd[,4],yminus=mean[,4]-sd[,4],bg='forestgreen', pch =22,
         type='b',add=T)
  errbar(time,mean[,5],yplus=mean[,5]+sd[,5],yminus=mean[,5]-sd[,5],bg='forestgreen', pch=23,
         type='b',add=T)
  errbar(time,mean[,6],yplus=mean[,6]+sd[,6],yminus=mean[,6]-sd[,6],bg='green', pch=24,
         type='b',add=T)
  legend('topleft',legend=c('NB4','NB6','NB7','NB7-22'), 
         pt.bg =c('forestgreen','forestgreen','forestgreen','green'),
         pch=c(21,22,23,24),bty='n')
  text(4,200,"NB4")
  
  
  # R
  errbar(time,mean[,7],yplus=mean[,7]+sd[,7],yminus=mean[,7]-sd[,7],
         log='y',type='b',ylim=c(119908.6,7787004.6), pch = 21, bg='navy',
         ylab="Cell Abundance (cells/mL)",
         xlab="Time (days)")
  errbar(time,mean[,9],yplus=mean[,9]+sd[,9],yminus=mean[,9]-sd[,9],bg='navy', pch =22,
         type='b',add=T)
  errbar(time,mean[,10],yplus=mean[,10]+sd[,10],yminus=mean[,10]-sd[,10],bg='navy', pch =23,
         type='b',add=T)
  errbar(time,mean[,11],yplus=mean[,11]+sd[,11],yminus=mean[,11]-sd[,11],bg='blue', pch=24,
         type='b',add=T)
  legend('topleft',legend=c('RA4','RA5','RD4','RD4-22'), 
         pt.bg =c('navy','navy','navy','blue'),
         pch=c(21,22,23,24),bty='n')
  text(4,200,"NB4")
  
  # R
  errbar(time,mean[,14],yplus=mean[,14]+sd[,14],yminus=mean[,14]-sd[,14],
         log='y',type='b',ylim=c(119908.6,7787004.6), pch =21, bg='firebrick',
         ylab="Cell Abundance (cells/mL)",
         xlab="Time (days)")
  errbar(time,mean[,15],yplus=mean[,15]+sd[,15],yminus=mean[,15]-sd[,15],bg='firebrick', pch = 22,
         type='b',add=T)
  legend('topleft',legend=c('YC5','YE5'), 
         pt.bg =c('firebrick','firebrick'),
         pch=c(21,22),bty='n')
  
  # R
  errbar(time,mean[,8],yplus=mean[,8]+sd[,8],yminus=mean[,8]-sd[,8],
         log='y',type='b',ylim=c(119908.6,7787004.6), pch =24,
         ylab="Cell Abundance (cells/mL)", bg= 'blue',
         xlab="Time (days)")
  errbar(time,mean[,12],yplus=mean[,12]+sd[,12],yminus=mean[,12]-sd[12],bg='green', pch = 24,
         type='b',add=T)
  errbar(time,mean[,13],yplus=mean[,13]+sd[,13],yminus=mean[,13]-sd[13],bg='orange', pch = 24,
         type='b',add=T)
  legend('topleft',legend=c('RD4 + RA5','RD4 + NB7', 'RD4 + YE5'), 
         pt.bg =c('blue','green','orange'),
         pch=24,bty='n')

  ## try with just average bacteria
read=read.csv(file="/Users/ahern/Desktop/bact_growth2.csv",header=T,row.names=1)
  
id=read$X
data=data.frame(read[,3:7])

  mean1=aggregate(data,by=list(id),FUN=mean,na.rm=T)
  mean=t(mean1[,2:6])
  sd1=aggregate(data,by=list(id),FUN=sd,na.rm=T)
  sd=t(sd1[2:6])
  # set time
  time=c(0,2,3,4,5)
  par(mfrow=c(1,2),mar=c(5,5,1,1))
  # Bact
errbar(time,mean[,1],yplus=mean[,1]+sd[,1],yminus=mean[,1]-sd[,1],
         log='y',type='b',ylim=c(103359.5,1682134.1),
         ylab="Cell Abundance (cells/mL)",
         xlab="Time (days)")
errbar(time,mean[,2],yplus=mean[,2]+sd[,2],yminus=mean[,2]-sd[,2],col='gray70',
         type='b',add=T)


read=read.csv(file="/Users/ahern/Desktop/bact_growth2.csv",header=T,row.names=1)

id=read$X
treatment=read$treat
data=data.frame(read[,3:7])

mean1=aggregate(data,by=list(treatment),FUN=mean,na.rm=T)
mean=t(mean1[,2:6])
sd1=aggregate(data,by=list(treatment),FUN=sd,na.rm=T)
sd=t(sd1[2:6])
# set time
time=c(0,2,3,4,5)
par(mfrow=c(1,2),mar=c(5,5,1,1))
# Bact
errbar(time,mean[,1],yplus=mean[,1]+sd[,1],yminus=mean[,1]-sd[,1],
       log='y',type='b',ylim=c(103359.5,1182134.1),
       ylab="Cell Abundance (cells/mL)",
       xlab="Time (days)")
errbar(time,mean[,2],yplus=mean[,2]+sd[,2],yminus=mean[,2]-sd[,2],col='gray70',
       type='b',add=T)



# NB
errbar(time,mean[,3],yplus=mean[,3]+sd[,3],yminus=mean[,3]-sd[,3],
       log='y',type='b',ylim=c(119908.6,7787004.6),
       ylab="Cell Abundance (cells/mL)",
       xlab="Time (days)")
errbar(time,mean[,4],yplus=mean[,4]+sd[,4],yminus=mean[,4]-sd[,4],col='gray70',
       type='b',add=T)
errbar(time,mean[,5],yplus=mean[,5]+sd[,5],yminus=mean[,5]-sd[,5],col='gray40',
       type='b',add=T)
errbar(time,mean[,6],yplus=mean[,6]+sd[,6],yminus=mean[,6]-sd[,6],col='gray40',
       type='b',add=T)
text(4,200,"NB4")


# R
errbar(time,mean[,7],yplus=mean[,7]+sd[,7],yminus=mean[,7]-sd[,7],
       log='y',type='b',ylim=c(119908.6,7787004.6),
       ylab="Cell Abundance (cells/mL)",
       xlab="Time (days)")
errbar(time,mean[,9],yplus=mean[,9]+sd[,9],yminus=mean[,9]-sd[,9],col='gray70',
       type='b',add=T)
errbar(time,mean[,10],yplus=mean[,10]+sd[,10],yminus=mean[,10]-sd[,10],col='gray40',
       type='b',add=T)
errbar(time,mean[,11],yplus=mean[,11]+sd[,11],yminus=mean[,11]-sd[,11],col='gray40',
       type='b',add=T)
text(4,200,"NB4")

# R
errbar(time,mean[,14],yplus=mean[,14]+sd[,14],yminus=mean[,14]-sd[,14],
       log='y',type='b',ylim=c(119908.6,7787004.6),
       ylab="Cell Abundance (cells/mL)",
       xlab="Time (days)")
errbar(time,mean[,15],yplus=mean[,15]+sd[,15],yminus=mean[,15]-sd[,15],col='gray70',
       type='b',add=T)



  
  
  
  

read=read.csv(file="/Users/ahern/Desktop/bact_gen.csv",header=T,row.names=1)
boxplot(read$gen_time1~read$treat,las=2)
treat=read$treat]
# bact diff?
bact=read$gen_time1[1:7]
a1=aov(read$gen_time1~read$x)
TukeyHSD(a1) # p= 0.0841159
boxplot(read$gen_time1~read$x)

# pops diff?
bact=read$gen_time1[8:37]
a1=aov(bact~pop[8:37])
summary(a1)
TukeyHSD(a1) # p= 0.5844931


# treat diff?
bact=read$gen_time1[8:37]
a1=aov(bact~treat[8:37])
summary(a1)
TukeyHSD(a1) # p= 0.5844931

boxplot(read$gen_2~read$treat)

# bact diff?
bact=read$gen_2[1:7]
a1=aov(read$gen_2~read$x)
TukeyHSD(a1) # p= 0.0841159
# pops diff?
bact=read$gen_2[8:37]
a1=aov(bact~pop[8:37])
TukeyHSD(a1) # p= 0.5844931
boxplot(bact~pop[8:37])
# treat diff?
bact=read$gen_2[8:37]
a1=aov(bact~treat[8:37])
TukeyHSD(a1) # p= 0.5844931
boxplot(bact~treat[8:37])





read=read.csv("/Users/ahern/Desktop/cg_certain_lefse.csv",header=T,row.names=1)
dim(read)
read=read.csv("/Users/ahern/Desktop/bact_by_lef.csv",header=T,row.names=1)


abund=((read/58270)*100)
colors=RColorBrewer::brewer.pal(8,'Spectral')

# par(mar=c(5,8,3,1),xpd=T,mfrow=c(2,1))
barplot(as.matrix(abund),ylim=c(0,40),names=rep("",81),space=0,border=NA,col=colors,
        ylab='Bacterial Class \nRelative Abundance',cex.lab=1,las=2,xlim=c(0,81))

barplot(as.matrix(abund[1:2,]),space=0,border=NA,col=colors,names=rep("",30),
        ylab='Bacterial Class \nRelative Abundance',cex.lab=1,las=2,xlim=c(0,30))



par(mar=c(20,10,10,1),xpd=T)
barplot(as.matrix(abund[1:3,]),las=2,space=0,border=NA,col=colors,
        ylab='Bacterial Family \nRelative Abundance',cex.lab=1.3,ylim=c(0,40))
box(which='plot')
segments(0,103,1.75,103,lwd=2)
text(.75,110,"T0-B",srt=0)
# T5s
segments(2,0,2,6,lwd=2,col='black')
segments(2.25,7,5.75,7,col='gray50',lwd=2)
text(4,110,"T5-B",srt=0)
# marches
segments(6,0,6,100,lwd=2,col='white')
segments(9,0,9,100,col='white')
segments(12,0,12,100,col='white')
segments(15,0,15,100,lwd=2,col='white')
segments(6.25,103,14.75,103,col='navy',lwd=2)
text(10.5,110,"CG1")
# mays
segments(18,0,18,100,col='white')
segments(21,0,21,100,lwd=2,col='white')
segments(15.25,103,20.75,103,col='firebrick',lwd=2)
text(18,110,"CG2")
# November
segments(24,0,24,100,lwd=1,col='white')
segments(27,0,27,100,col='white')
segments(30,0,30,100,lwd=2,col='white')
segments(21.25,103,30,103,col='forestgreen',lwd=2)
text(25.5,110,"CG3")



par(mar=c(20,10,10,1),xpd=T)
barplot(as.matrix(abund[3:8,]),las=2,space=0,border=NA,col=colors,
        ylab='Bacterial Family \nRelative Abundance',cex.lab=1.3,ylim=c(0,40))
box(which='plot')
segments(0,103,1.75,103,lwd=2)
text(.75,110,"T0-B",srt=0)
# T5s
segments(2,0,2,6,lwd=2,col='black')
segments(2.25,7,5.75,7,col='gray50',lwd=2)
text(4,110,"T5-B",srt=0)
# marches
segments(6,0,6,100,lwd=2,col='white')
segments(9,0,9,100,col='white')
segments(12,0,12,100,col='white')
segments(15,0,15,100,lwd=2,col='white')
segments(6.25,103,14.75,103,col='navy',lwd=2)
text(10.5,110,"CG1")
# mays
segments(18,0,18,100,col='white')
segments(21,0,21,100,lwd=2,col='white')
segments(15.25,103,20.75,103,col='firebrick',lwd=2)
text(18,110,"CG2")
# November
segments(24,0,24,100,lwd=1,col='white')
segments(27,0,27,100,col='white')
segments(30,0,30,100,lwd=2,col='white')
segments(21.25,103,30,103,col='forestgreen',lwd=2)
text(25.5,110,"CG3")


colors=RColorBrewer::brewer.pal(12,'Set3')

par(mar=c(5,5,5,1),mfrow=c(2,2),xpd=T)
barplot(as.matrix(abund[1:3,]),las=2,space=0,border=NA,col=colors,
        ylab='Bacterial OTU Relative Abundance',cex.lab=1.3,ylim=c(0,10))
{box(which='plot')
segments(0,101,1.75,101,lwd=2)
text(.75,105,"T0-B",srt=0,cex=2)
# T5s
segments(2,0,2,100,lwd=2,col='black')
segments(2.25,101,5.75,101,col='gray50',lwd=2)
text(4,105,"T5-B",srt=0,cex=2)
# marches
segments(6,0,6,100,lwd=2,col='black')
segments(9,0,9,100,col='gray50')
segments(12,0,12,100,col='gray50')
segments(15,0,15,100,lwd=2,col='gray50')
segments(6.25,101,14.75,101,col='navy',lwd=2)
text(10.5,105,"CG1",cex=2)
# mays
segments(18,0,18,100,col='gray50')
segments(21,0,21,100,lwd=2,col='black')
segments(15.25,101,20.75,101,col='firebrick',lwd=2)
text(18,105,"CG2",cex=2)
# November
segments(24,0,24,100,lwd=1,col='gray50')
segments(27,0,27,100,col='gray50')
segments(30,0,30,100,lwd=2,col='gray50')
segments(21.25,101,30,101,col='forestgreen',lwd=2)
text(25.5,105,"CG3",cex=2)}
# Bact
barplot(as.matrix(abund[4:6,]),las=2,space=0,border=NA,col=colors[4:6],
        ylab='Bacterial Family \nRelative Abundance',cex.lab=1.3,ylim=c(0,10))
{box(which='plot')
  segments(0,101,1.75,101,lwd=2)
  text(.75,105,"T0-B",srt=0,cex=2)
  # T5s
  segments(2,0,2,100,lwd=2,col='black')
  segments(2.25,101,5.75,101,col='gray50',lwd=2)
  text(4,105,"T5-B",srt=0,cex=2)
  # marches
  segments(6,0,6,100,lwd=2,col='black')
  segments(9,0,9,100,col='gray50')
  segments(12,0,12,100,col='gray50')
  segments(15,0,15,100,lwd=2,col='gray50')
  segments(6.25,101,14.75,101,col='navy',lwd=2)
  text(10.5,105,"CG1",cex=2)
  # mays
  segments(18,0,18,100,col='gray50')
  segments(21,0,21,100,lwd=2,col='black')
  segments(15.25,101,20.75,101,col='firebrick',lwd=2)
  text(18,105,"CG2",cex=2)
  # November
  segments(24,0,24,100,lwd=1,col='gray50')
  segments(27,0,27,100,col='gray50')
  segments(30,0,30,100,lwd=2,col='gray50')
  segments(21.25,101,30,101,col='forestgreen',lwd=2)
  text(25.5,105,"CG3",cex=2)}
barplot(as.matrix(abund[7:9,]),las=2,space=0,border=NA,col=colors[7:9],
        ylab='Bacterial Family \nRelative Abundance',cex.lab=1.3,ylim=c(0,10))
{box(which='plot')
  segments(0,101,1.75,101,lwd=2)
  text(.75,105,"T0-B",srt=0,cex=2)
  # T5s
  segments(2,0,2,100,lwd=2,col='black')
  segments(2.25,101,5.75,101,col='gray50',lwd=2)
  text(4,105,"T5-B",srt=0,cex=2)
  # marches
  segments(6,0,6,100,lwd=2,col='black')
  segments(9,0,9,100,col='gray50')
  segments(12,0,12,100,col='gray50')
  segments(15,0,15,100,lwd=2,col='gray50')
  segments(6.25,101,14.75,101,col='navy',lwd=2)
  text(10.5,105,"CG1",cex=2)
  # mays
  segments(18,0,18,100,col='gray50')
  segments(21,0,21,100,lwd=2,col='black')
  segments(15.25,101,20.75,101,col='firebrick',lwd=2)
  text(18,105,"CG2",cex=2)
  # November
  segments(24,0,24,100,lwd=1,col='gray50')
  segments(27,0,27,100,col='gray50')
  segments(30,0,30,100,lwd=2,col='gray50')
  segments(21.25,101,30,101,col='forestgreen',lwd=2)
  text(25.5,105,"CG3",cex=2)}
barplot(as.matrix(abund[13:15,]),las=2,space=0,border=NA,col=colors[10:12],
        ylab='Bacterial Family \nRelative Abundance',cex.lab=1.3,ylim=c(0,10))
{box(which='plot')
  segments(0,101,1.75,101,lwd=2)
  text(.75,105,"T0-B",srt=0,cex=2)
  # T5s
  segments(2,0,2,100,lwd=2,col='black')
  segments(2.25,101,5.75,101,col='gray50',lwd=2)
  text(4,105,"T5-B",srt=0,cex=2)
  # marches
  segments(6,0,6,100,lwd=2,col='black')
  segments(9,0,9,100,col='gray50')
  segments(12,0,12,100,col='gray50')
  segments(15,0,15,100,lwd=2,col='gray50')
  segments(6.25,101,14.75,101,col='navy',lwd=2)
  text(10.5,105,"CG1",cex=2)
  # mays
  segments(18,0,18,100,col='gray50')
  segments(21,0,21,100,lwd=2,col='black')
  segments(15.25,101,20.75,101,col='firebrick',lwd=2)
  text(18,105,"CG2",cex=2)
  # November
  segments(24,0,24,100,lwd=1,col='gray50')
  segments(27,0,27,100,col='gray50')
  segments(30,0,30,100,lwd=2,col='gray50')
  segments(21.25,101,30,101,col='forestgreen',lwd=2)
  text(25.5,105,"CG3",cex=2)}



```

# 16S data

```
## Warning: package 'DESeq2' was built under R version 4.3.3
```

```
## Warning: package 'GenomeInfoDb' was built under R version 4.3.3
```

```
## [1] 1581   94
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1581 taxa and 94 samples ]
## sample_data() Sample Data:       [ 94 samples by 29 sample variables ]
## tax_table()   Taxonomy Table:    [ 1581 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1581 tips and 1580 internal nodes ]
```

```
## [1] TRUE
```


# General Alpha statistics 

```r
cg_filt=subset_samples(cg, keep_2026 =='yes')
s=specnumber(t(otu_table(cg_filt)))

# spec number for 
filter=subset_samples(cg_filt,clonal2=='monoclonal' & definition =='assemblage')
filter
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1574 taxa and 12 samples ]
## sample_data() Sample Data:       [ 12 samples by 29 sample variables ]
## tax_table()   Taxonomy Table:    [ 1574 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1574 tips and 1573 internal nodes ]
```

```r
s2=specnumber(t(otu_table(filter)))
mean(s2)
```

```
## [1] 277.3333
```

```r
sd(s2)
```

```
## [1] 45.88787
```

```r
filter=subset_samples(cg_filt,clonal2=='multiclonal' & definition =='assemblage')
s=specnumber(t(otu_table(filter)))
mean(s)
```

```
## [1] 309.5556
```

```r
sd(s)
```

```
## [1] 30.02962
```

```r
t.test(s, s2)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  s and s2
## t = 1.9407, df = 18.745, p-value = 0.06749
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -2.56108 67.00552
## sample estimates:
## mean of x mean of y 
##  309.5556  277.3333
```

# Beta diversity and ordinations

## variance stabilizing transformation 


```r
library(CoDaSeq)
library(DESeq2)
cg_filt <- filter_taxa(cg_filt, function(x) sum(x) > 1, prune = TRUE)

########################################################################
# variance stabilizing transformation 
# from 16S stats 
ddst=phyloseq_to_deseq2(cg_filt,~Treatment3)
# 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
trial.deseq=estimateSizeFactors(ddst,geoMeans = apply(counts(ddst), 1, gm_mean))
data.vst <- DESeq2::varianceStabilizingTransformation(trial.deseq, blind = T)
# Apply variance stabilizing transformation
# Extract transformed OTU table
fieldEco.vstMat <- assay(data.vst)
# And put into my phloseq object
OTU = otu_table(fieldEco.vstMat, taxa_are_rows=T)
comp_vst=phyloseq(OTU,map,tax2,tree)
comp_clr=comp_vst
```

## all

### pca

```r
library(CoDaSeq)
library(DESeq2)
cg_filt2=subset_samples(cg_filt)
########################################################################
# variance stabilizing transformation 
# from 16S stats 
ddst=phyloseq_to_deseq2(cg_filt2,~Treatment3)
# 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
trial.deseq=estimateSizeFactors(ddst,geoMeans = apply(counts(ddst), 1, gm_mean))
data.vst <- DESeq2::varianceStabilizingTransformation(trial.deseq, blind = T)
# Apply variance stabilizing transformation
# Extract transformed OTU table
fieldEco.vstMat <- assay(data.vst)
# And put into my phloseq object
OTU = otu_table(fieldEco.vstMat, taxa_are_rows=T)
comp_vst=phyloseq(OTU,map,tax2,tree)
comp_clr=comp_vst
p=prcomp(t(otu_table(comp_clr)))

{plot(p$x[,1],p$x[,2],pch=sample_data(comp_clr)$pch,
     bg=as.character(sample_data(comp_clr)$col),
     xlab = "PCA 43%", ylab= "PCA 11%")
ordiellipse(p$x,
     #    group=as.factor(sample_data(comp_clr)$Treatment),
          group=as.factor(sample_data(comp_clr)$Pop3),
           col=c('black', 'navy', 'firebrick','forestgreen'),

     #     col=c('navy', 'firebrick', 'forestgreen','black'),
           kind ='sd',conf=0.8,
           label=TRUE)

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/vst-2-1.png)<!-- -->

### pcoa

```r
otu2=otu_table(comp_clr)
euc=phyloseq::distance(t(otu2),'euclidean')
p=prcomp(euc,scale=T,center=T)

{
  
  par(mar=c(10,
            5,1,1),xpd=FALSE)
  plot(p$x[,1],p$x[,2],pch=sample_data(comp_clr)$pch,
     bg=as.character(sample_data(comp_clr)$col2),
     xlab = "PCoA1 68%", 
     ylab= "PCoA2 21%",
     yaxt='n')
  axis(2, at=c(-2,0,2))
ordiellipse(p$x,
     #    group=as.factor(sample_data(comp_clr)$Treatment),
          group=as.factor(sample_data(comp_clr)$Pop3),
           col=c('black', 'navy', 'firebrick','forestgreen'),

     #     col=c('navy', 'firebrick', 'forestgreen','black'),
           kind ='sd',conf=0.8,
           label=TRUE)
legend(-2, -5,
       legend=c("RA5","RD4", "YC5", "NB7", "RD4 + RA5", "RD4 + YE5", "RD4+NB7"),
       bty='n',
  #    col=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),
       pt.bg=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),

       pch=c(22,23,22,23,21,22,23))

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/vst-2-pcoa-1.png)<!-- -->

```r
anosim(euc, grouping=as.factor(sample_data(comp_clr)$Pop3))
```

```
## 
## Call:
## anosim(x = euc, grouping = as.factor(sample_data(comp_clr)$Pop3)) 
## Dissimilarity: euclidean 
## 
## ANOSIM statistic R: 0.1382 
##       Significance: 0.01 
## 
## Permutation: free
## Number of permutations: 999
```

## microbiomes only 

### pca

```r
library(CoDaSeq)
library(DESeq2)
cg_filt2=subset_samples(cg_filt, extract=="WGA" & definition =="microbiome")
########################################################################
# variance stabilizing transformation 
# from 16S stats 
ddst=phyloseq_to_deseq2(cg_filt2,~Treatment3)
# 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
trial.deseq=estimateSizeFactors(ddst,geoMeans = apply(counts(ddst), 1, gm_mean))
data.vst <- DESeq2::varianceStabilizingTransformation(trial.deseq, blind = T)
# Apply variance stabilizing transformation
# Extract transformed OTU table
fieldEco.vstMat <- assay(data.vst)
# And put into my phloseq object
OTU = otu_table(fieldEco.vstMat, taxa_are_rows=T)
comp_vst=phyloseq(OTU,map,tax2,tree)
comp_clr=comp_vst
p=prcomp(t(otu_table(comp_clr)))
# summary(p)

{plot(p$x[,1],p$x[,2],pch=sample_data(comp_clr)$pch,
     bg=as.character(sample_data(comp_clr)$col),
     xlab = "PCA 19 %", ylab= "PCA 7 %")
ordiellipse(p$x,
     #    group=as.factor(sample_data(comp_clr)$Treatment),
          group=as.factor(sample_data(comp_clr)$Pop3),
           col=c('black', 'navy', 'firebrick','forestgreen'),

     #     col=c('navy', 'firebrick', 'forestgreen','black'),
           kind ='sd',conf=0.8,
           label=TRUE)

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/vst-3-1.png)<!-- -->

### pcoa


```r
otu2=otu_table(comp_clr)
euc=phyloseq::distance(t(otu2),'euclidean')
p=prcomp(euc,scale=T,center=T)

# anosim(euc, grouping=as.factor(sample_data(comp_clr)$Pop3))
# anosim(euc, grouping=as.factor(sample_data(comp_clr)$diatom_control))

{
  
  par(mar=c(15,
            5,1,1),xpd=FALSE)
  plot(p$x[,1],p$x[,2],pch=sample_data(comp_clr)$pch,
     bg=as.character(sample_data(comp_clr)$col2),
     xlab = "PCoA1 36 %",
     ylab= "PCoA2 10 %", 
     yaxt='n', xaxt='n')
axis(2, at=c(-2,0,4))
axis(1, at=c(-4,0,5,10))
ordiellipse(p$x,
     #    group=as.factor(sample_data(comp_clr)$Treatment),
          group=as.factor(sample_data(comp_clr)$Pop3),
           col=c('black', 'navy', 'firebrick','forestgreen'),

     #     col=c('navy', 'firebrick', 'forestgreen','black'),
           kind ='sd',conf=0.8,
           label=TRUE)
legend(-2, -5,
       legend=c("RA5","RD4", "YC5", "NB7", "RD4 + RA5", "RD4 + YE5", "RD4+NB7"),
       bty='n',
  #    col=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),
       pt.bg=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),

       pch=c(22,23,22,23,21,22,23))

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/vst-4-1.png)<!-- -->

```r
anosim(euc, grouping=as.factor(sample_data(comp_clr)$Pop3))
```

```
## 
## Call:
## anosim(x = euc, grouping = as.factor(sample_data(comp_clr)$Pop3)) 
## Dissimilarity: euclidean 
## 
## ANOSIM statistic R: 0.2457 
##       Significance: 0.001 
## 
## Permutation: free
## Number of permutations: 999
```

```r
# cg_filta=subset_samples(comp_clr, clonal=='monoclonal')
# otu2=otu_table(cg_filta)
# euc=phyloseq::distance(t(otu2),'euclidean')
# anosim(euc, grouping=as.factor(sample_data(cg_filta)$Pop3))
# anosim(euc, grouping=as.factor(sample_data(cg_filta)$Pop3))
# adonis2(euc~as.factor(sample_data(cg_filta)$Pop3))
```


## assemblages 

### pca

```r
library(CoDaSeq)
library(DESeq2)
cg_filt2=subset_samples(cg_filt, definition =='assemblage')
########################################################################
# variance stabilizing transformation 
# from 16S stats 
ddst=phyloseq_to_deseq2(cg_filt2,~Treatment3)
# 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
trial.deseq=estimateSizeFactors(ddst,geoMeans = apply(counts(ddst), 1, gm_mean))
data.vst <- DESeq2::varianceStabilizingTransformation(trial.deseq, blind = T)
# Apply variance stabilizing transformation
# Extract transformed OTU table
fieldEco.vstMat <- assay(data.vst)
# And put into my phloseq object
OTU = otu_table(fieldEco.vstMat, taxa_are_rows=T)
comp_vst=phyloseq(OTU,map,tax2,tree)
comp_clr=comp_vst


p=prcomp(t(otu_table(comp_clr)))
# summary(p)

{plot(p$x[,1],p$x[,2],pch=sample_data(comp_clr)$pch,
     bg=as.character(sample_data(comp_clr)$col),
     xlab = "PCA 28 %", ylab= "PCA 11 %")
ordiellipse(p$x,
     #    group=as.factor(sample_data(comp_clr)$Treatment),
          group=as.factor(sample_data(comp_clr)$Pop3),
           col=c('black', 'navy', 'firebrick','forestgreen'),

     #     col=c('navy', 'firebrick', 'forestgreen','black'),
           kind ='sd',conf=0.8,
           label=TRUE)

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/vst-5-1.png)<!-- -->

### pcoa

```r
otu2=otu_table(comp_clr)
euc=phyloseq::distance(t(otu2),'euclidean')
p=prcomp(euc,scale=T,center=T)
# summary(p)

{
  
  par(mar=c(10,
            5,1,1),xpd=FALSE)
  plot(p$x[,1],p$x[,2],pch=sample_data(comp_clr)$pch,
     bg=as.character(sample_data(comp_clr)$col2),
     xlab = "PCoA1 41 %",
     ylab= "PCoA2 10 %",
     yaxt='n')
  axis(2, at=c(-2,0,2))
ordiellipse(p$x,
     #    group=as.factor(sample_data(comp_clr)$Treatment),
          group=as.factor(sample_data(comp_clr)$Pop3),
           col=c('black', 'navy', 'firebrick','forestgreen'),

     #     col=c('navy', 'firebrick', 'forestgreen','black'),
           kind ='sd',conf=0.8,
           label=TRUE)
legend(-2, -5,
       legend=c("RA5","RD4", "YC5", "NB7", "RD4 + RA5", "RD4 + YE5", "RD4+NB7"),
       bty='n',
  #    col=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),
       pt.bg=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),

       pch=c(22,23,22,23,21,22,23))

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/vst-6-1.png)<!-- -->

## assemblages with control 

```r
library(CoDaSeq)
library(DESeq2)
cg_filt2=subset_samples(cg_filt, definition =='assemblage' | definition == "inoculum_T5")
########################################################################
# variance stabilizing transformation 
# from 16S stats 
ddst=phyloseq_to_deseq2(cg_filt2,~Treatment3)
# 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
trial.deseq=estimateSizeFactors(ddst,geoMeans = apply(counts(ddst), 1, gm_mean))
data.vst <- DESeq2::varianceStabilizingTransformation(trial.deseq, blind = T)
# Apply variance stabilizing transformation
# Extract transformed OTU table
fieldEco.vstMat <- assay(data.vst)
# And put into my phloseq object
OTU = otu_table(fieldEco.vstMat, taxa_are_rows=T)
comp_vst=phyloseq(OTU,map,tax2,tree)
comp_clr=comp_vst

otu2=otu_table(comp_clr)
euc=phyloseq::distance(t(otu2),'euclidean')
p=prcomp(euc,scale=T,center=T)
# summary(p)


# anosim(euc, grouping=as.factor(sample_data(comp_clr)$Pop3), strata=sample_data(comp_clr)$Treatment2)
# anosim(euc, grouping=(sample_data(comp_clr)$diatom_control))
# adonis2(euc ~ (sample_data(comp_clr)$diatom_control))

{
  
  par(mar=c(15,
            5,1,1),xpd=FALSE)
  plot(p$x[,1],p$x[,2],pch=sample_data(comp_clr)$pch,
     bg=as.character(sample_data(comp_clr)$col2),
     xlab = "PCoA1 41.11 %",
     ylab= "PCoA2 10.18 %",
     yaxt='n')
  axis(2, at=c(-2,0,2))
ordiellipse(p$x,
     #    group=as.factor(sample_data(comp_clr)$Treatment),
          group=as.factor(sample_data(comp_clr)$Pop3),
           col=c('black', 'navy', 'firebrick','forestgreen'),

     #     col=c('navy', 'firebrick', 'forestgreen','black'),
           kind ='sd',conf=0.8,
           label=TRUE)
legend(-2, -5,
       legend=c("RA5","RD4", "YC5", "NB7", "RD4 + RA5", "RD4 + YE5", "RD4+NB7"),
       bty='n',
  #    col=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),
       pt.bg=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),

       pch=c(22,23,22,23,21,22,23))

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/vst-7-1.png)<!-- -->

# philr 

### all


```r
library(philr)
sub=subset_samples(cg_filt, extract=="WGA" & definition =="microbiome")

GP <- transform_sample_counts(sub, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "strain_Otu2116/Kingdom_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
tree <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')

gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(GP, 'PCoA', distance=gp.dist)
gp.p=prcomp(gp.dist)

{
  par(mar=c(15,
            5,1,1),xpd=F)
  plot(gp.p$x[,1],gp.p$x[,2],pch=sample_data(sub)$pch,
     bg=as.character(sample_data(sub)$col2),
       xlab = "PCoA1 38.52 %",
     ylab= "PCoA2 19.90 %")
# axis(2, at=c(-2,0,2))
ordiellipse(gp.p$x,
    #    group=as.factor(sample_data(comp_clr)$Treatment),
         group=as.factor(sample_data(sub)$Pop3),
           col=c('black', 'navy', 'firebrick','forestgreen'),

     #     col=c('navy', 'firebrick', 'forestgreen','black'),
           kind ='sd',conf=0.8,
           label=TRUE)
legend(-20,-40,
        legend=c("RA5","RD4", "YC5", "NB7", "RD4 + RA5", "RD4 + YE5", "RD4 + NB7"),
       bty='n',
  #    col=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),
       pt.bg=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),

       pch=c(22,23,22,23,21,22,23))

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/philr-1-1.png)<!-- -->

### microbiomes monoclonal

```r
library(philr)
sub=subset_samples(cg_filt, extract=="WGA" & definition =="microbiome" & clonal =='monoclonal')

GP <- transform_sample_counts(sub, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "strain_Otu2116/Kingdom_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
tree <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')

gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(GP, 'PCoA', distance=gp.dist)
gp.p=prcomp(gp.dist)

{
  par(mar=c(15,
            5,1,1),xpd=F)
  plot(gp.p$x[,1],gp.p$x[,2],pch=sample_data(sub)$pch,
     bg=as.character(sample_data(sub)$col2),
       xlab = "PCoA1 39.68 %",
     ylab= "PCoA2 19.46 %")
# axis(2, at=c(-2,0,2))
ordiellipse(gp.p$x,
    #    group=as.factor(sample_data(comp_clr)$Treatment),
         group=as.factor(sample_data(sub)$Pop3),
           col=c('black', 'navy', 'firebrick','forestgreen'),

     #     col=c('navy', 'firebrick', 'forestgreen','black'),
           kind ='sd',conf=0.8,
           label=TRUE)
legend(-20,-40,
        legend=c("RA5","RD4", "YC5", "NB7", "RD4 + RA5", "RD4 + YE5", "RD4 + NB7"),
       bty='n',
  #    col=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),
       pt.bg=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),

       pch=c(22,23,22,23,21,22,23))

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/philr-2-1.png)<!-- -->


### assemblage control 

```r
library(philr)
sub=subset_samples(cg_filt, definition =='assemblage' | definition == "inoculum_T5")
GP <- transform_sample_counts(sub, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "strain_Otu2116/Kingdom_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
tree <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')

gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(GP, 'PCoA', distance=gp.dist)
gp.p=prcomp(gp.dist)

{plot(gp.pcoa$vectors[,1],gp.pcoa$vectors[,2],pch=sample_data(sub)$pch,
     bg=as.character(sample_data(sub)$col),
     xlab = "PCoA 1 53 %", ylab= "PCoA 2 28%")
ordiellipse(gp.pcoa$vectors,
           group=as.factor(sample_data(sub)$Treatment),label=T,
           col=c('gray70','gray50','forestgreen','green','navy','blue'),
           kind ='sd',conf=0.8)}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/philr-3-1.png)<!-- -->

### assemblage only 

```r
library(philr)
sub=subset_samples(cg_filt, definition =='assemblage')
GP <- transform_sample_counts(sub, function(x) x+1)
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "strain_Otu2116/Kingdom_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
tree <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')

gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(GP, 'PCoA', distance=gp.dist)
gp.p=prcomp(gp.dist)

{
  par(mar=c(15,
            5,1,1),xpd=F)
  plot(gp.p$x[,1],gp.p$x[,2],pch=sample_data(sub)$pch,
     bg=as.character(sample_data(sub)$col2),
       xlab = "PCoA1 74  %",
     ylab= "PCoA2 8 %")
# axis(2, at=c(-2,0,2))
ordiellipse(gp.p$x,
    #    group=as.factor(sample_data(comp_clr)$Treatment),
         group=as.factor(sample_data(sub)$Pop3),
           col=c('black', 'navy', 'firebrick','forestgreen'),

     #     col=c('navy', 'firebrick', 'forestgreen','black'),
           kind ='sd',conf=0.8,
           label=TRUE)
legend(-20,-40,
        legend=c("RA5","RD4", "YC5", "NB7", "RD4 + RA5", "RD4 + YE5", "RD4 + NB7"),
       bty='n',
  #    col=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),
       pt.bg=c("navy",'navy', 'firebrick', 'forestgreen', 'black','black','black'),

       pch=c(22,23,22,23,21,22,23))

}
```

![](/Users/oliviaahern/Documents/GitHub/Diatom_Microbiome/docs/index_files/figure-html/philr-4-1.png)<!-- -->
