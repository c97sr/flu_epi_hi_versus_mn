rm(list=ls(all=TRUE))
# setwd("~/Dropbox/shares/KMW_SR/project2a")

# Third party libraries
require("xtable")
require("mgcv")
require("fitdistrplus")

# Package source for SRileyIDD is available as tar.gz from http://goo.gl/O9Qdai
# Then "R CMD install SRileyIDD_0.1.tar.gz" in any dir at the command prompt
# Or in R install.packages("SRileyIDD_0.1.tar.gz",repos=NULL,type="source")
# Or git clone https://github.com/c97sr/SRileyIDD.git SRileyIDD 
# then R CMD install SRileyIDD
require("SRileyIDD")

##############

# Code for to setup the data tables

##############

# Load the full dataset
dat_full <- read.csv("./hkiss_main_v3.csv")

# Add variables that will be needed in all subsets
dat_full$Ag1 <- cut(dat_full$Age_1,breaks=c(0,19,44,65,110),labels=c("0_18","19_44","45_64","65_110"))
dat_full$yMN <- dat_full$PH1N1_com_result_tbase
dat_full$yHI <- (dat_full$H1N1.T2 / dat_full$H1N1.T1) > 2.1
dat_full$t1HIlog <- log2(dat_full$H1N1.T1 / 5)
dat_full$t2HIlog <- log2(dat_full$H1N1.T2 / 5)
dat_full$logboost <- dat_full$t2HIlog - dat_full$t1HIlog
dat_full$lbcl <- ifelse(dat_full$logboost < 3,NA,dat_full$logboost)
dat_full$lbcr <- ifelse(dat_full$logboost < 3,2,dat_full$logboost)
dat_full$mn1 <- match(dat_full$ph1n1.T1.Raw,c("01:10","01:20","01:40"))
dat_full$mn2 <- match(dat_full$ph1n1.T2.Raw,c("01:10","01:20","01:40"))
dat_full$mn1known <- !is.na(dat_full$mn1)
dat_full$mn2known <- !is.na(dat_full$mn2)
dat_full$mn12known <- dat_full$mn1known & dat_full$mn2known

# Cut to the 770 with good data for this and the fields we will use 
fields_selected <- c(
  "New_id","H1N1.T1","H1N1.T2","H1N1MN.T1","H1N1MN.T2",
  "PH1N1_com_result_tbase","Age_1","District.5","child_hh",
  "Ag1","yMN","yHI","t1HIlog","t2HIlog","logboost","lbcr","lbcl",
  "ph1n1.T1.Raw","ph1n1.T2.Raw","mn1","mn2","mn1known","mn2known",
  "mn12known"
) 
dat_large <- dat_full[	
    (!is.na(dat_full$PH1N1_com_result_tbase)),  
    fields_selected
]
dat <- dat_full[	
  (!is.na(dat_full$PH1N1_com_result_tbase)) & 
  !is.na(dat_full$H1N1.T1), 
  fields_selected
]

##############

# Code for results paragraph 1 about the raw comparison

##############

# Including one for a discrepancy
dat$disc <- ((dat$yMN) & (!dat$yHI)) | ((!dat$yMN) & (dat$yHI)) 

# R1 Make a table of the discrepancies 
dat$Tab1Type <- rep("NULL",dim(dat)[1])
dat$Tab1Type[dat$yMN & dat$yHI] <- "1 MN+ HI+"
dat$Tab1Type[(!dat$yMN) & (!dat$yHI)] <- "2 MN- HI- "
dat$Tab1Type[(dat$yMN) & (!dat$yHI)] <- "3 MN+ HI-"
dat$Tab1Type[(!dat$yMN) & (dat$yHI)] <- "4 MN- HI+"
xtable(addmargins(table(dat$Ag1,dat$Tab1Type), margin=c(1,2)),align="rccccc",display=c("d","d","d","d","d","d"))
dim(dat)

# Test the logistic regression model for discrepanciy by age
# No difference here at the moment
mod_disc_smooth <- gam(disc ~ s(Age_1), family=binomial, data=dat)
mod_disc_cat <- gam(disc ~ as.factor(Ag1), family=binomial, data=dat)
summary(mod_disc_smooth)
summary(mod_disc_cat)

##############

# Code for results for paragraph 2 about the different models
# Micro neutralization seems to give a much better explination

##############

# Start with the final model from PLoS Med
mod_LINEAR_FULL_MN <- gam(
  yMN ~ as.numeric(Age_1) + as.factor(child_hh),
  data=dat_full,family=binomial
)
mod_LINEAR_SMALL_MN <- gam(
  yMN ~ as.numeric(Age_1) + as.factor(child_hh),
  data=dat,family=binomial
)
mod_LINEAR_SMALL_HI <- gam(
  yHI ~ as.numeric(Age_1) + as.factor(child_hh),
  data=dat,family=binomial
)
mod_SMOOTH_SMALL_MN <- gam(
  yMN ~ s(Age_1) + as.factor(child_hh),
  data=dat,family=binomial
)
mod_SMOOTH_SMALL_HI <- gam(
  yHI ~ s(Age_1) + as.factor(child_hh),
  data=dat,family=binomial
)

# Make a file called "modelsum.csv" that can be used as the basis for 
# an excel file for the latex table
table2 <- summarise.glm(
  list(
    mod_LINEAR_FULL_MN,mod_LINEAR_SMALL_MN,mod_LINEAR_SMALL_HI,
    mod_SMOOTH_SMALL_MN,mod_SMOOTH_SMALL_HI),
  transpose=FALSE,
  file="./table2_raw.csv"
)

# A plot to show that the smooth mdoel for MN is just the linear model
# Probbaly not included in the thesis but maybe in the 
plot(mod_SMOOTH_SMALL_MN)
dev.off()

pdf.figure.proof(file="./fig1_proof.pdf")
plot(mod_SMOOTH_SMALL_HI)
dev.off()

##############

# Code for results paragraph 3 about boosting

##############

# All measures for which both HI and MN are known
indexMaskaa <- c(which(dat$mn1known),which(dat$mn2known))
dat_aa <- data.frame(
  id=dat$New_id[indexMaskaa],
  age=dat$Age_1[indexMaskaa],
  ag=dat$Ag1[indexMaskaa],
  hi=c(dat$t1HIlog[dat$mn1known],dat$t2HIlog[dat$mn2known]),
  mn=c(dat$mn1[dat$mn1known],dat$mn2[dat$mn2known]) 
)
dim(dat_aa)

pdf.figure.proof(file="./fig2_proof.pdf")
plot(jitter(dat_aa$mn),jitter(dat_aa$hi))
dev.off()

sum(dat$mn1known)
sum(dat$mn2known)
sum(dat$mn2known) + sum(dat$mn1known)
addmargins(table(dat_aa$hi,dat_aa$mn))

aa_mod1 <- gam(hi ~ as.numeric(mn),data=dat_aa,family=poisson)
aa_mod2 <- gam(hi ~ as.numeric(mn) + as.factor(ag),data=dat_aa,family=poisson)
aa_mod3 <- gam(hi ~ as.numeric(mn),data=dat_aa[dat$mn1known,],family=poisson)
aa_mod4 <- gam(hi ~ as.numeric(mn) + as.factor(ag),data=dat_aa[dat$mn1known,],family=poisson)
aa_mod5 <- gam(hi ~ as.numeric(mn),data=dat_aa[dat$mn2known,],family=poisson)
aa_mod6 <- gam(hi ~ as.numeric(mn) + as.factor(ag),data=dat_aa[dat$mn2known,],family=poisson)

table3 <- summarise.glm(
  list(aa_mod1,aa_mod2,aa_mod3,aa_mod4,aa_mod5,aa_mod6),
  transpose=FALSE,
  file="./table3_raw.csv"
)
