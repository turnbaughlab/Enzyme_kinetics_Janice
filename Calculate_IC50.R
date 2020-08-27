## =========================================================================
## Calculating IC50 values 
## Authors: Janice Goh
## =========================================================================

## Settings 
setwd("/Users/jgoh/Box Sync/For_Huy/")
library(dplyr)
library(MESS)
library(drc)
library(reshape)

## load data 
data <- read.csv("data.format.csv")

## calculate mean, SD and CV for each compound
## do note i am doing this a little old school so please forgive me 
data.melt <- melt(data, id.vars = c("Chem_ID", "Concentration..uM."))
data.stats <- data.melt %>% group_by(Chem_ID, Concentration..uM.) %>%
              summarise(mean = mean(value),
                        SD = sd(value),
                        CV = sd(value)/mean(value))
## Get IC50 values
data.split <- split(data.stats, data.stats$Chem_ID)
IC50.all <- NULL
pdf(file = "IC50.graphs.pdf")
for(chem.idx in names(data.split)){
  df <- data.split[[chem.idx]]
  IC50.plot <- drm(df$mean~df$Concentration..uM., data = df, fct = LL.4 (fixed = c(NA, NA, 100, NA))) 
  ## By defining the third parameter, I have set the max possible response to 100. 
  ## Set 100 to NA if you think otherwise
  ## if you expect a full DRC i.e. you want compounds to go from 100 to 0, then set the 2nd parameter to 0.
  ## if you do this, some curves will fail to converge because there is no actual inhibition 
  ## extract IC50 value
  IC50.val <- IC50.plot$fit$par[3]
  IC50.all <- cbind(IC50.all, IC50.val)
  ## plot IC50 curves 
  plot(IC50.plot, 
      xlab = ("Drug concentration (uM)"), 
      ylab= ("Response (%)"),
      ylim = c(0, 150),
      main = chem.idx)
  arrows(df$Concentration..uM., df$mean-df$SD,
         df$Concentration..uM., df$mean+df$SD, code=3, length=0.02, angle = 90)
  mtext(side=3, line= 0, text = paste0("IC50=", round(IC50.val, digits = 3)))
}
dev.off()
colnames(IC50.all) <- names(data.split)

## save IC50 results as csv file
write.csv(IC50.all, file = "IC50.values.csv")
