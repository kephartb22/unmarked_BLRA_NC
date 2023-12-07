# unmarked_BLRA_NC
single-season species occupancy modelling for Black Rails in NC


##1.) Loading packages

library(unmarked)
library(MuMIn)
library(dplyr)
library(readr)
library(AICcmodavg)

##2.) Loading in data

data <- read_csv("occupancy model 120623/model_covariates_formatted_v6.csv")
View(data)


##3.) Establish observation covariates

lunar <- data[ ,c("lunar.1.scaled", "lunar.2.scaled", "lunar.3.scaled", "lunar.4.scaled", "lunar.5.scaled")]

julian.day <- data[ ,c("julianday.1.scaled", "julianday.2.scaled", "julianday.3.scaled", "julianday.4.scaled", "julianday.5.scaled")]

wind <- data[ ,c("wind.1.scaled", "wind.2.scaled", "wind.3.scaled", "wind.4.scaled", "wind.5.scaled")]


##4.) Establish site covariates

sitecovariates_df <- data[ ,c("juncus.centered", "spart.centered", "water.centered", "tree.centered", "iva.centered", "phrag.centered", "VO_standard", "DS_standard", "burn_scaled", "elev.scaled")]
siteCovs <- data.frame(sitecovariates_df)


##5.) Establish y

y <- data[,4:8]

##5.) Create model frame

umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = list(julian.day = julian.day, lunar = lunar, wind = wind))
summary(umf)


##6.) Summarize detection histories

detHist(umf)


##7.) Create models

###null model

summary(null<-occu(~1~1,data=umf))



backTransform(null, 'state')


###global model

summary(full <- occu(formula = ~julian.day + lunar + wind          ~juncus.centered + spart.centered + water.centered + tree.centered + burn_scaled + VO_standard + DS_standard + elev.scaled, data=umf))



modelList_full <- dredge(full,
                    rank = "AICc")

model_df_full <- as.data.frame(modelList_full)


###site only

summary(full_SITE <- occu(formula = ~1          ~juncus.centered + spart.centered + water.centered + tree.centered + burn_scaled + VO_standard + DS_standard + elev.scaled, data=umf, control = list(maxit = 500)))

####control = list(maxit =) used because model does not converge without increased maximum iterations

modelList_site <- dredge(full_SITE,
                    rank = "AICc")

model_df_SITE <- as.data.frame(modelList_site)


###obs only

summary(full_OBS <- occu(formula = ~julian.day + lunar + wind          ~1, data=umf, control = list(maxit = 500)))

modelList_obs <- dredge(full_OBS,
                    rank = "AICc")

model_df_OBS <- as.data.frame(modelList_obs)



###model averaging top models with delta AICc <2

###call in top models

mod_rank1 <- occu(formula = ~1          ~tree.centered + water.centered + elev.scaled, data=umf)
mod_rank2 <- occu(formula = ~wind          ~tree.centered + water.centered + elev.scaled, data=umf)
mod_rank3 <- occu(formula = ~1          ~tree.centered + water.centered, data=umf)
mod_rank4 <- occu(formula = ~1          ~tree.centered + water.centered + elev.scaled + spart.centered, data=umf)
mod_rank5 <- occu(formula = ~1          ~elev.scaled + juncus.centered + tree.centered + water.centered, data=umf)
mod_rank6 <- occu(formula = ~1          ~burn_scaled + elev.scaled + tree.centered + water.centered, data=umf)



model_top_list <- list(mod_rank1, mod_rank2, mod_rank3, mod_rank4, mod_rank5, mod_rank6)

model_avg <- model.avg(model_top_list, rank = AICc)
summary(model_avg)




###MB GOF for top model

occ_gof <- mb.gof.test(mod_rank1, nsim= 1000, plot.hist = FALSE)
occ_gof$chisq.table <- NULL
print(occ_gof)


