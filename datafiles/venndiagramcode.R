# ven code 
## venn diagrams


### nb7
```{r, echo=FALSE}

sample_data(cg)$Treatment

nb7=subset_samples(cg, Treatment=="NB7" | Treatment=="NB7_Blank")
sample_data(nb7)

agg=tax_glom(nb7, taxrank="species")
library(grid)
library(VennDiagram)

spec1=specnumber(t(otu_table(subset_samples(agg, Treatment_Multi=="NB7_Filter"))))
nb7 = merge_samples((agg), "Treatment_Multi")
Bacterioplankton <- prune_samples("NB7_Blank_WGA", nb7)
Bacterioplankton <- filter_taxa(Bacterioplankton, function(x) sum(x) > 1, prune = TRUE)
Bacterioplankton=colnames((otu_table(Bacterioplankton)))
Bacterioplankton
Assemblage <- prune_samples("NB7_Filter", nb7)
Assemblage <- filter_taxa(Assemblage, function(x) sum(x) > 1, prune = TRUE)
Assemblage=colnames((otu_table(Assemblage)))
Assemblage
Microbiome <- prune_samples("NB7_WGA", nb7)
Microbiome <- filter_taxa(Microbiome, function(x) sum(x) > 1, prune = TRUE)
Microbiome=colnames((otu_table(Microbiome)))
Microbiome


calculate.overlap(
  x=list("Bacterioplankton"=Bacterioplankton, 
         "Microbiome"=Microbiome, 
         "Assemblage"=Assemblage))

shared=intersect(row.names(Bacterioplankton),row.names(Microbiome))
length(shared)
nb7_bacterio=Bacterioplankton
nb7_micro=Microbiome
nb7_assem=Assemblage

venn1=venn.diagram(x=list(nb7_bacterio,
                          nb7_micro,
                          nb7_assem),
                   category.names=c("Bacterioplankton",
                                    "Microbiome"
                                    ,"Assemblage"),
                   filename = NULL,
                   main="PopC NB7",
                   #  main.cex=1.4,
                   height=1,width=1,
                   cat.cex=1,
                   scaled = TRUE,
                   fontfamily=rep("sans", 7),
                   cat.fontfamily=rep("sans", 3))

cowplot::plot_grid(venn1)

ven=venn.diagram(x=list(nb7_bacterio,
                        nb7_micro,
                        nb7_assem),
                 category.names=c("Bacterioplankton",
                                  "Microbiome"
                                  ,"Assemblage"),
                 filename = NULL,
                 main="PopC NB7",
                 #  main.cex=1.4,
                 height=1,width=1,
                 cat.cex=1,
                 scaled = TRUE,
                 fontfamily=rep("sans", 7),
                 cat.fontfamily=rep("sans", 3))

cowplot::plot_grid(ven)


cam=calculate.overlap(x=list(nb7_bacterio,
                             nb7_micro,
                             nb7_assem))

names(cam) <- c("all",
                "Bacterio-Micro",
                "Bacterio-Assemb",
                "Micro-Assemb",
                "Bacterio", "Micro", "Assemb")


cam
tax_mat <- as(tax_table(nb7), "matrix")
taxa_to_keep=c("Otu650")
subset_tax_mat <- tax_mat[rownames(tax_mat) %in% taxa_to_keep, ]
subset_tax_mat

taxa_to_keep=c("Otu37")
subset_tax_mat <- tax_mat[rownames(tax_mat) %in% taxa_to_keep, ]
subset_tax_mat

cowplot::plot_grid(venn1)

length(cam$all)
length(cam$Bacterio)
length(cam$Micro)
length(cam$Assemb)
length(cam$`Bacterio-Micro`)
length(cam$`Bacterio-Assemb`)
length(cam$`Micro-Assemb`)


```

### RA5
```{r, echo=FALSE}
sample_data(cg)$Treatment

RA5=subset_samples(cg, Treatment=="RA5" | Treatment=="RA5_Blank")
sample_data(RA5)

agg=tax_glom(RA5, taxrank="species")
library(grid)
library(VennDiagram)

RA5 = merge_samples((agg), "Treatment_Multi")
Bacterioplankton <- prune_samples("RA5_Blank_WGA", RA5)
Bacterioplankton <- filter_taxa(Bacterioplankton, function(x) sum(x) > 1, prune = TRUE)
Bacterioplankton=colnames((otu_table(Bacterioplankton)))
Bacterioplankton
Assemblage <- prune_samples("RA5_Filter", RA5)
Assemblage <- filter_taxa(Assemblage, function(x) sum(x) > 1, prune = TRUE)
Assemblage=colnames((otu_table(Assemblage)))
Assemblage
Microbiome <- prune_samples("RA5_WGA", RA5)
Microbiome <- filter_taxa(Microbiome, function(x) sum(x) > 1, prune = TRUE)
Microbiome=colnames((otu_table(Microbiome)))
Microbiome

venn_ra5=venn.diagram(x=list(Bacterioplankton,Microbiome, Assemblage),
                      category.names=c("Bacterioplankton","Microbiome","Assemblage"),
                      filename = NULL,
                      main="PopA RA5",
                      main.cex=2,
                      main.fontface	='bold')

ra5_bacterio=Bacterioplankton
ra5_micro=Microbiome
ra5_assem=Assemblage


venn_ra5=venn.diagram(x=list(ra5_bacterio,ra5_micro, ra5_assem),
                      category.names=c("Bacterioplankton","Microbiome","Assemblage"),
                      filename = NULL,
                      main="PopA RA5",
                      main.cex=1.4,
                      height=1,width=1,
                      cat.cex=1,
                      scaled = TRUE,
                      fontfamily=rep("sans", 7),
                      cat.fontfamily=rep("sans", 3),
                      ext.text = TRUE,
                      ext.line.lwd = 2,
                      ext.dist = -0.15,
                      ext.length = 0.9,
                      ext.pos = -4)

cowplot::plot_grid(venn1, venn_ra5)


cam=calculate.overlap(x=list(ra5_bacterio,
                             ra5_micro,
                             ra5_assem))

names(cam) <- c("all",
                "Bacterio-Micro",
                "Bacterio-Assemb",
                "Micro-Assemb",
                "Bacterio", "Micro", "Assemb")
cam
tax_mat <- as(tax_table(RA5), "matrix")
taxa_to_keep=cam$Bacterio
subset_tax_mat <- tax_mat[rownames(tax_mat) %in% taxa_to_keep, ]
subset_tax_mat
taxa_to_keep=cam$Micro
subset_tax_mat <- tax_mat[rownames(tax_mat) %in% taxa_to_keep, ]
subset_tax_mat



length(cam$all)
length(cam$Bacterio)
length(cam$Micro)
length(cam$Assemb)
length(cam$`Bacterio-Micro`)
length(cam$`Bacterio-Assemb`)
length(cam$`Micro-Assemb`)
```

### RD4

```{r, echo=FALSE}

sample_data(cg)$Treatment

RD4=subset_samples(cg, Treatment=="RD4" | Treatment=="RD4_Blank")
sample_data(RD4)

agg=tax_glom(RD4, taxrank="species")
library(grid)
library(VennDiagram)

RD4 = merge_samples((agg), "Treatment_Multi")
Bacterioplankton <- prune_samples("RD4_Blank_WGA", RD4)
Bacterioplankton <- filter_taxa(Bacterioplankton, function(x) sum(x) > 1, prune = TRUE)
Bacterioplankton=colnames((otu_table(Bacterioplankton)))
Bacterioplankton
Assemblage <- prune_samples("RD4_Filter", RD4)
Assemblage <- filter_taxa(Assemblage, function(x) sum(x) > 1, prune = TRUE)
Assemblage=colnames((otu_table(Assemblage)))
Assemblage
Microbiome <- prune_samples("RD4_WGA", RD4)
Microbiome <- filter_taxa(Microbiome, function(x) sum(x) > 1, prune = TRUE)
Microbiome=colnames((otu_table(Microbiome)))
Microbiome


rd4_bacterio=Bacterioplankton
rd4_micro=Microbiome
rd4_assem=Assemblage

venn_rd4=venn.diagram(x=list(rd4_bacterio,rd4_micro, rd4_assem),
                      category.names=c("Bacterioplankton","Microbiome","Assemblage"),
                      filename = NULL,
                      main="PopA RD4",
                      main.cex=1.4,
                      fontfamily=rep("sans", 7),
                      cat.fontfamily=rep("sans", 3),
                      height=1,width=1,
                      cat.cex=1,
                      scaled = TRUE,
                      ext.text = TRUE,
                      ext.line.lwd = 2,
                      ext.dist = -0.15,
                      ext.length = 0.9,
                      ext.pos = -4)
cowplot::plot_grid(venn_rd4)

cowplot::plot_grid(venn1, venn_ra5,venn_rd4,scale=0.8)


cam=calculate.overlap(x=list(rd4_bacterio,
                             rd4_micro,
                             rd4_assem))

names(cam) <- c("all",
                "Bacterio-Micro",
                "Bacterio-Assemb",
                "Micro-Assemb",
                "Bacterio", "Micro", "Assemb")
cam
tax_mat <- as(tax_table(RD4), "matrix")
taxa_to_keep=cam$Bacterio
subset_tax_mat <- tax_mat[rownames(tax_mat) %in% taxa_to_keep, ]
subset_tax_mat


taxa_to_keep=cam$Micro
subset_tax_mat <- tax_mat[rownames(tax_mat) %in% taxa_to_keep, ]
subset_tax_mat



length(cam$all)
length(cam$Bacterio)
length(cam$Micro)
length(cam$Assemb)
length(cam$`Bacterio-Micro`)
length(cam$`Bacterio-Assemb`)
length(cam$`Micro-Assemb`)

```


### YE5

```{r, echo=FALSE}

#sample_data(cg)$Treatment

YE5=subset_samples(cg, Treatment=="YE5" | Treatment=="YE5_Blank")
#sample_data(YE5)

agg=tax_glom(YE5, taxrank="species")
library(grid)
library(VennDiagram)

YE5 = merge_samples((agg), "Treatment_Multi")
sample_data(YE5)
Bacterioplankton <- prune_samples("YE5_Blank_WGA", YE5)
Bacterioplankton <- filter_taxa(Bacterioplankton, function(x) sum(x) > 1, prune = TRUE)
Bacterioplankton=colnames((otu_table(Bacterioplankton)))
Bacterioplankton
Assemblage <- prune_samples("YE5_Filter", YE5)
Assemblage <- filter_taxa(Assemblage, function(x) sum(x) > 1, prune = TRUE)
Assemblage=colnames((otu_table(Assemblage)))
Assemblage
Microbiome <- prune_samples("YE5_WGA", YE5)
Microbiome <- filter_taxa(Microbiome, function(x) sum(x) > 1, prune = TRUE)
Microbiome=colnames((otu_table(Microbiome)))
Microbiome

venn_ye5=venn.diagram(x=list(Bacterioplankton,Microbiome, Assemblage),
                      category.names=c("Bacterioplankton","Microbiome","Assemblage"),
                      filename = NULL,
                      main="PopB YE5",
                      fontfamily=rep("sans", 7),
                      cat.fontfamily=rep("sans", 3),
                      main.cex=1.4,
                      height=1,width=1,
                      cat.cex=1,
                      scaled = TRUE,
                      ext.text = TRUE,
                      ext.line.lwd = 2,
                      ext.dist = -0.15,
                      ext.length = 0.9,
                      ext.pos = -4)
cowplot::plot_grid(venn_ye5)

cowplot::plot_grid(venn1, venn_ra5,venn_rd4,venn_ye5,
                   scale=0.8)
cowplot::plot_grid(venn_ra5,venn_rd4,venn_ye5,venn1,
                   nrow=1,
                   scale=0.5)



cam=calculate.overlap(x=list(Bacterioplankton,
                             Microbiome,
                             Assemblage))

names(cam) <- c("all",
                "Bacterio-Micro",
                "Bacterio-Assemb",
                "Micro-Assemb",
                "Bacterio", "Micro", "Assemb")
cam
tax_mat <- as(tax_table(YE5), "matrix")
taxa_to_keep=c("Otu752")
subset_tax_mat <- tax_mat[rownames(tax_mat) %in% taxa_to_keep, ]
subset_tax_mat


taxa_to_keep=c("Otu276","Otu650","Otu1301", "Otu738","Otu827")
subset_tax_mat <- tax_mat[rownames(tax_mat) %in% taxa_to_keep, ]
subset_tax_mat



length(cam$all)
length(cam$Bacterio)
length(cam$Micro)
length(cam$Assemb)
length(cam$`Bacterio-Micro`)
length(cam$`Bacterio-Assemb`)
length(cam$`Micro-Assemb`)
```





