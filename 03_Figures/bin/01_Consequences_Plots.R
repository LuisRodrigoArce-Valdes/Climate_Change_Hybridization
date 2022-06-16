#!/usr/bin/Rscript
# 01_Consequences_Plots.R
# Made by Luis Rodrigo Arce Valdés, to plot genetic distance for each evolutionary consequence
rm(list = ls())

# Calling libraries
library(dplyr)
library(ggplot2)
library(tidyr)
# library(scales)

# Reading input files ####
COIs <- read.delim("../../02_COIs/results/01_Genetic_Distances.tsv", header = T)
Ugly.Data <- read.delim("../../01_Review/results/Database_Full.txt", header = T)
Tidy.Data <- read.delim("../../01_Review/results/TidyConsequences.tsv", header = T)

# First we will include genetic distances into our tables, and fixing bees subspecies
gsub(" ","_", paste0(Ugly.Data$Sp1,"_X_",Ugly.Data$Sp2)) -> Ugly.Data$Cross
Ugly.Data[Ugly.Data=="Apis_mellifera_scutellata_X_Apis_mellifera_ssp"] <- "Apis_mellifera_scutellata_X_Apis_mellifera_ligustica"

gsub(" ","_", paste0(Tidy.Data$Sp1,"_X_",Tidy.Data$Sp2)) -> Tidy.Data$Cross
Tidy.Data[Tidy.Data=="Apis_mellifera_scutellata_X_Apis_mellifera_ssp"] <- "Apis_mellifera_scutellata_X_Apis_mellifera_ligustica"

# Joining
left_join(Ugly.Data, COIs, by="Cross") %>% select(!Cross) -> Ugly.Data
write.table(Ugly.Data, "../results/01_Full_Table.tsv", sep = "\t", quote = F, row.names = F)
rm(Ugly.Data)

# Now for the tidy data
left_join(Tidy.Data, COIs, by="Cross") %>% select(!Cross) -> Tidy.Data

# Violin plot ####
# First plot; Violin plots (removing species that still have not hybridized):
violins <- Tidy.Data[Tidy.Data$Type!="Expected Range Expansion And Hybridisation",]
violins <- violins[complete.cases(violins),]

# Categoryzing outcomes
factor(violins$Outcome, levels = unique(violins$Outcome),
       labels = c("Neutral","Positive","Negative",
                  "Positive","Negative","Positive",
                  "Neutral","Negative","Positive")) -> violins$Effect
violins$Effect <- factor(violins$Effect, levels = c("Positive","Neutral","Negative"))

# Sorting
violins <- violins[order(violins$Effect, violins$Outcome),]
violins[violins=="Species Fussion"] <- "Species Fusion"

# factoring outcomes
violins$Outcome <- factor(violins$Outcome, levels = rev(unique(violins$Outcome)))

# Plotting
png("../results/02_Violins.png", width = 32, height = 16, units = "cm", res = 300)
ggplot(violins) +
  geom_violin(aes(x=Distance, y=Outcome, color=Effect), fill="grey80", scale = "width", size=1) +
  geom_point(aes(x=Distance, y=Outcome, fill=Order), size=2.5, shape=21) +
  scale_color_manual(values = c('#1b9e77','#7570b3','#d95f02')) +
  scale_fill_manual(values = c('#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6')) +
  stat_summary(aes(x=Distance, y=Outcome), fun=function(y) median(y), geom='point', color="white", size=1, shape=3) +
  theme_classic() +
  labs(x = "Genetic Distance") +
  theme(text = element_text(size = 20, family = "serif"),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 13),
        legend.title = element_blank())
dev.off()

# Statistical testing:
sink("../results/03_Violins_testing.txt", append = F, split = T)
kruskal.test(Distance ~ Outcome, data = violins)
for(i in c("none","bonferroni","BY")) {
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(pairwise.wilcox.test(violins$Distance, violins$Outcome, p.adjust.method = i))
}
sink()

# Radarchart code ####
radarchart <- function(df, axistype=0, seg=4, pty=16, pcol=1:8, plty=1:6, plwd=1,
                       pdensity=NULL, pangle=45, pfcol=NA, cglty=3, cglwd=1,
                       cglcol="navy", axislabcol="blue", title="", maxmin=TRUE,
                       na.itp=TRUE, centerzero=FALSE, vlabels=NULL, vlcex=NULL,
                       caxislabels=NULL, calcex=NULL,
                       paxislabels=NULL, palcex=NULL, ...) {
  if (!is.data.frame(df)) { cat("The data must be given as dataframe.\n"); return() }
  if ((n <- length(df))<3) { cat("The number of variables must be 3 or more.\n"); return() }
  if (maxmin==FALSE) { # when the dataframe does not include max and min as the top 2 rows.
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type="n", frame.plot=FALSE, axes=FALSE, 
       xlab="", ylab="", main=title, asp=1, ...) # define x-y coordinates without any plot
  theta <- seq(90, 450, length=n+1)*pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) { # complementary guide lines, dotted navy line by default
    polygon(xx*(i+CGap)/(seg+CGap), yy*(i+CGap)/(seg+CGap), lty=cglty, lwd=cglwd, border=cglcol)
    if (axistype==1|axistype==3) CAXISLABELS <- paste(i/seg*100,"(%)")
    if (axistype==4|axistype==5) CAXISLABELS <- sprintf("%3.2f",i/seg)
    if (!is.null(caxislabels)&(i<length(caxislabels))) CAXISLABELS <- caxislabels[i+1]
    if (axistype==1|axistype==3|axistype==4|axistype==5) {
      if (is.null(calcex)) text(-0.05, (i+CGap)/(seg+CGap), CAXISLABELS, col=axislabcol) else
        text(-0.05, (i+CGap)/(seg+CGap), CAXISLABELS, col=axislabcol, cex=calcex)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)
  }
  else {
    arrows(xx/(seg+CGap), yy/(seg+CGap), xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)
  }
  PAXISLABELS <- df[1,1:n]
  if (!is.null(paxislabels)) PAXISLABELS <- paxislabels
  if (axistype==2|axistype==3|axistype==5) {
    if (is.null(palcex)) text(xx[1:n], yy[1:n], PAXISLABELS, col=axislabcol) else
      text(xx[1:n], yy[1:n], PAXISLABELS, col=axislabcol, cex=palcex)
  }
  VLABELS <- colnames(df)
  
  if (!is.null(vlabels)) VLABELS <- vlabels
  
  
  
  
  
  ##--------------------------------------------------
  ## Modified by Killbill-(Me)
  ##--------------------------------------------------
  # Main code:
  # if (is.null(vlcex)) text(xx*1.2, yy*1.2, VLABELS) else
  #   text(xx*1.2, yy*1.2, VLABELS, cex=vlcex, adj=adjVec)
  
  
  # Modified code:
  # Create a variable that round 'xx' value to 0 and 1 for non zero and 0.5 for 0 values.
  adjVec <- ifelse(round(xx) < 0, 1, ifelse(round(xx) > 0, 0, 0.5))
  
  #apply 'adjVec' variable to "adj" parameters of text.
  
  for (i in seq_along(xx)){
    if (is.null(vlcex)) text(xx[i]*1.1, yy[i]*1.1, VLABELS[i], adj=adjVec[i]) else
      text(xx[i]*1.1, yy[i]*1.1, VLABELS[i], cex=vlcex, adj=adjVec[i])
  }
  
  ##-------------------------------------------------
  ## End
  ##-------------------------------------------------
  
  
  
  
  series <- length(df[[1]])
  SX <- series-2
  if (length(pty) < SX) { ptys <- rep(pty, SX) } else { ptys <- pty }
  if (length(pcol) < SX) { pcols <- rep(pcol, SX) } else { pcols <- pcol }
  if (length(plty) < SX) { pltys <- rep(plty, SX) } else { pltys <- plty }
  if (length(plwd) < SX) { plwds <- rep(plwd, SX) } else { plwds <- plwd }
  if (length(pdensity) < SX) { pdensities <- rep(pdensity, SX) } else { pdensities <- pdensity }
  if (length(pangle) < SX) { pangles <- rep(pangle, SX)} else { pangles <- pangle }
  if (length(pfcol) < SX) { pfcols <- rep(pfcol, SX) } else { pfcols <- pfcol }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg+CGap)+(df[i,]-df[2,])/(df[1,]-df[2,])*seg/(seg+CGap)
    if (sum(!is.na(df[i,]))<3) { cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n",i,df[i,])) # for too many NA's (1.2.2012)
    } else {
      for (j in 1:n) {
        if (is.na(df[i, j])) { # how to treat NA
          if (na.itp) { # treat NA using interpolation
            left <- ifelse(j>1, j-1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left>1, left-1, n)
            }
            right <- ifelse(j<n, j+1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right<n, right+1, 1)
            }
            xxleft <- xx[left]*CGap/(seg+CGap)+xx[left]*(df[i,left]-df[2,left])/(df[1,left]-df[2,left])*seg/(seg+CGap)
            yyleft <- yy[left]*CGap/(seg+CGap)+yy[left]*(df[i,left]-df[2,left])/(df[1,left]-df[2,left])*seg/(seg+CGap)
            xxright <- xx[right]*CGap/(seg+CGap)+xx[right]*(df[i,right]-df[2,right])/(df[1,right]-df[2,right])*seg/(seg+CGap)
            yyright <- yy[right]*CGap/(seg+CGap)+yy[right]*(df[i,right]-df[2,right])/(df[1,right]-df[2,right])*seg/(seg+CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft; yytmp <- yyleft;
              xxleft <- xxright; yyleft <- yyright;
              xxright <- xxtmp; yyright <- yytmp;
            }
            xxs[j] <- xx[j]*(yyleft*xxright-yyright*xxleft)/(yy[j]*(xxright-xxleft)-xx[j]*(yyright-yyleft))
            yys[j] <- (yy[j]/xx[j])*xxs[j]
          } else { # treat NA as zero (origin)
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j]*CGap/(seg+CGap)+xx[j]*(df[i, j]-df[2, j])/(df[1, j]-df[2, j])*seg/(seg+CGap)
          yys[j] <- yy[j]*CGap/(seg+CGap)+yy[j]*(df[i, j]-df[2, j])/(df[1, j]-df[2, j])*seg/(seg+CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty=pltys[i-2], lwd=plwds[i-2], border=pcols[i-2], col=pfcols[i-2])
      } else {
        polygon(xxs, yys, lty=pltys[i-2], lwd=plwds[i-2], border=pcols[i-2], 
                density=pdensities[i-2], angle=pangles[i-2], col=pfcols[i-2])
      }
      points(xx*scale, yy*scale, pch=ptys[i-2], col=pcols[i-2])
    }
  }
}
# Spider Plots ####

# $$$$$$$$$ To edit categories names $$$$$$$ #
Tidy.Data[Tidy.Data=="Climate Change Induced Hybridisation"] <- "Climate Change Induced Hybridization"
Tidy.Data[Tidy.Data=="Hybrid Zone Dynamics Dependant On Climate Change"] <- "Hybrid Zone Dynamics\nDependent On Climate Change"
Tidy.Data[Tidy.Data=="Increased Expected Hybridisation Due To Climate Change"] <- "Hybrid Zone Dynamics\nDependent On Climate Change"
Tidy.Data[Tidy.Data=="Expected Range Expansion And Hybridisation"] <- "Expected Range Expansion\nAnd Hybridization"
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #

# Creating a dataframe per plot
spiders <- list()
spiders[["Outcomes"]] <- as.data.frame(table(Tidy.Data$Outcome))
Tidy.Data %>% select(!Outcome) %>% 
  unique() -> Tidy.Data
spiders[["Order"]] <- as.data.frame(table(Tidy.Data$Order))
spiders[["Type"]] <- as.data.frame(table(Tidy.Data$Type))

# Now for Regions we will need to tidy a little bit
data.frame(Regions = Tidy.Data$Region) %>% 
  separate(Regions, into = c("A","B","C"), sep = ", ", fill = "right") %>% 
  gather("Extra","Region",1:3) %>% 
  select(Region) -> Regions

# Removing NAs
Regions <- Regions[complete.cases(Regions),]
spiders[["Regions"]] <- as.data.frame(table(Regions))

# Preparing dataframes
for(i in 1:length(spiders)){
  row.names(spiders[[i]]) <- spiders[[i]][,1]
  spiders[[i]][,1] <- 10
  spiders[[i]] <- as.data.frame.matrix(t(spiders[[i]]))
  spiders[[i]] <- spiders[[i]][-1,]
  spiders[[i]] <- rbind(rep(max(spiders[[i]][1,]), ncol(spiders[[i]])), rep(0, ncol(spiders[[i]])), spiders[[i]])
  spiders[[i]] <- spiders[[i]][,c(1,ncol(spiders[[i]]):2)]
}

# Editing names of outcomes orders
spiders$Outcomes <- spiders$Outcomes[,c(7,2,5,9,4,3,6,8,1)]
colnames(spiders$Outcomes)[5] <- "Rare\nHybridization"

# Changing order
spiders <- spiders[c("Type","Outcomes", "Order", "Regions")]

# Color vector
colors_border=c(rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9), rgb(0.7,0.5,0.1,0.9), rgb(0.4,0.5,0.3,0.9) )
colors_in=c(rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4), rgb(0.7,0.5,0.1,0.4), rgb(0.4,0.5,0.3,0.4) )

# Spider plot
for(i in 1:length(spiders)) {
  png(paste0("../results/04_Spider_Plot_",i,".png"), width = 45.5, height = 26, units = "cm", res = 300)
  par(family="serif")
  radarchart(spiders[[i]], axistype = 1,
             # custom polygon
             pcol=colors_border[i] , pfcol=colors_in[i] , plwd=2 , plty=1, caxislabels=seq(0,spiders[[i]][1,1], length.out=5), calcex=1.5,
             #custom the grid
             cglcol="black", cglty=1, axislabcol="black", cglwd=0.8,
             #custom labels
             vlcex=2.5)
  dev.off()
}

