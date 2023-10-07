install.packages("tableone")
library(tidyverse)
library(tableone)
library(survival)
library(broom)
stroke <- read.csv("stroke_rtw.csv")

# Date variables of date ------------------------------------------------
stroke$randomisation_dte<-as.Date(stroke$randomisation_dte)
stroke$rtw_dte<-as.Date(stroke$rtw_dte)
stroke$exit_assessment_dte<-as.Date(stroke$exit_assessment_dte)

# clean database stroke ---------------------------------------------
# clean database: step1: calculate interval free of outcome 
stroke <- stroke %>%
  # calculate rtw_date - randomisation_date (unit=year)
  mutate(rtw_random = ifelse(rtw_flg==1, rtw_dte - randomisation_dte, exit_assessment_dte-randomisation_dte)) %>%
  mutate(rtw_random = round(rtw_random/365.25, 2))
  
# clean database: step2: imputing missing interval free of outcome  
# calculate the avg_dte of rtw: rtw_dte_avg_ESSVR and rtw_dte_avg_USUAL
rtw_rand_avg_ESSVR <- stroke %>% filter(rtw_flg==1 & alloc=="ESSVR")
rtw_rand_avg_ESSVR <-  mean(rtw_rand_avg_ESSVR$rtw_random, na.rm=TRUE)

rtw_rand_avg_USUAL <- stroke %>% filter(rtw_flg==1 & alloc=="Usual Care")
rtw_rand_avg_USUAL <- mean(rtw_rand_avg_USUAL$rtw_random, na.rm=TRUE)

stroke <- stroke %>%
  # impute missing rtw_random for rtw_flg==1 & alloc=="ESSVR"
  mutate(rtw_random=if_else(is.na(rtw_random) & rtw_flg==1 & alloc=="ESSVR",rtw_rand_avg_ESSVR ,rtw_random)) %>%
  # impute missing rtw_random for rtw_flg==1 & alloc=="Usual Care"
  mutate(rtw_random=if_else(is.na(rtw_random) & rtw_flg==1 & alloc=="Usual Care",rtw_rand_avg_USUAL, rtw_random)) 

# factor categorical variables ------------------------------------------
stroke$essvr_complete_flg <- factor(stroke$essvr_complete_flg)
stroke$rtw_flg <- factor(stroke$rtw_flg)
stroke$alloc <- factor(stroke$alloc, levels=c("Usual Care", "ESSVR"))

# create time split dataset for timesplit cox regression model
tcuts <- c(0.5)
stroke_splite <-survSplit(Surv(rtw_random, rtw_flg)~., data=stroke, cut = tcuts, episode="ftime_group")
stroke_splite$ftime_group <- factor(stroke_splite$ftime_group)


# descriptive analysis ------------------------------------------------
ESSVR <- stroke %>% filter(!is.na(essvr_complete_flg))
attri_names <- dput(names(stroke))
attri_descri <- c("sex", "age", "region", "work_status_pre", "hpw_pre", "stroke_severity", "rtw_flg", "health_score")

# descriptive analysis 1 for ESSVR group, stratified by essvr_complete --
Descri_ana_ESSVR <- CreateTableOne(data=ESSVR, vars=attri_descri, strata="essvr_complete_flg")

summary(Descri_ana_ESSVR) # check if there is missing data

Descri_table_ESSVR <- print(Descri_ana_ESSVR, showAllLevels = TRUE) # show all levels of categorical variables

write.csv(Descri_table_ESSVR, file="Descri_table_ESSVR.csv")

# descriptive analysis 2 for all participants, stratified by alloc ------
Descri_ana_stroke <- CreateTableOne(data=stroke, vars=attri_descri, strata="alloc")

summary(Descri_ana_stroke) # check if there is missing data

Descri_table_stroke <- print(Descri_ana_stroke, showAllLevels = TRUE) # show all levels of categorical variables

write.csv(Descri_table_stroke, file="Descri_table_stroke.csv")

# --------------------------------------------------------------------
# AIM 1: Effect of work_status_pre, hpw_pre and stroke_severity on completing the ESSVR

# simple logistic regression model ------------------------------------
fit1_work_status <- glm(essvr_complete_flg ~ work_status_pre, family = binomial, data = ESSVR)
fit1_hpw <- glm(essvr_complete_flg ~ hpw_pre, family = binomial, data = ESSVR)
fit1_stroke_sev <- glm(essvr_complete_flg ~ stroke_severity, family = binomial, data = ESSVR)

# simple logistic regression: fit1_work_status on ESSVR completion
exp(coef(fit1_work_status))
exp(confint(fit1_work_status))
summary(fit1_work_status)$coefficients

# simple logistic regression: fit1_hpw on ESSVR completion
exp(coef(fit1_hpw))
exp(confint(fit1_hpw))
summary(fit1_hpw)$coefficients

# simple logistic regression: fit1_stroke on ESSVR completion
exp(coef(fit1_stroke_sev))
exp(confint(fit1_stroke_sev))
summary(fit1_stroke_sev)$coefficients

# multiple logistic regression model adjusted for age and sex --------
# work_status_pre and stroke_severity are only 2 significant variables
fit2_work_status <- glm(essvr_complete_flg ~ work_status_pre + sex + age, family = binomial, data = ESSVR)
fit2_hpw <- glm(essvr_complete_flg ~ hpw_pre + sex + age, family = binomial, data = ESSVR)
fit2_stroke_sev <- glm(essvr_complete_flg ~ stroke_severity + sex + age, family = binomial, data = ESSVR)

# multiple logistic regression: fit2_work_status on ESSVR completion
exp(coef(fit2_work_status))
exp(confint(fit2_work_status))
summary(fit2_work_status)$coefficients

# multiple logistic regression: fit2_hpw on ESSVR completion
exp(coef(fit2_hpw))
exp(confint(fit2_hpw))
summary(fit2_hpw)$coefficients

# multiple logistic regression: fit2_stroke on ESSVR completion
exp(coef(fit2_stroke_sev))
exp(confint(fit2_stroke_sev))
summary(fit2_stroke_sev)$coefficients

# --------------------------------------------------------------------
# AIM 2-1: Whether ESSVR affected return to work post-stroke(sub-group by sex) -----------------------------------------------------------------
# Effect of ESSVR on returning to work - follow up for 1 year 
fit3_rtw_ESSVR_alloc <- coxph(Surv(rtw_random, as.numeric(rtw_flg)) ~ alloc, data=stroke)

# examine the proportional hazards assumption
install.packages("survminer")
library(survminer)
library(survival)
# cumulative hazard plot of returning to work 
survest <- survfit(Surv(rtw_random, as.numeric(rtw_flg))~alloc, data=stroke)

ggsurvplot(survest, fun = "cumhaz", conf.int = TRUE, data=stroke, xlab="year", title="Cumulative hazard of returning to work")

# log-rank test

survdiff(formula=Surv(rtw_random, as.numeric(rtw_flg))~alloc, data=stroke)

# cox regression -------------------
survminer::ggcoxzph(cox.zph(fit3_rtw_ESSVR_alloc))


summary(fit3_rtw_ESSVR_alloc)
# summary(fit3_rtw_ESSVR_alloc_male)
# summary(fit3_rtw_ESSVR_alloc_female)

# AIM 2-2: Effect modification 
fit3_rtw_ESSVR_alloc_sex <- coxph(Surv(rtw_random, as.numeric(rtw_flg)) ~ alloc*sex, data=stroke)
summary(fit3_rtw_ESSVR_alloc_sex)

fit3_rtw_ESSVR_alloc_age <- coxph(Surv(rtw_random, as.numeric(rtw_flg)) ~ alloc*age, data=stroke)
summary(fit3_rtw_ESSVR_alloc_age)

fit3_rtw_ESSVR_alloc_work_status_pre <- coxph(Surv(rtw_random, as.numeric(rtw_flg)) ~ alloc*work_status_pre, data=stroke)
summary(fit3_rtw_ESSVR_alloc_work_status_pre)

fit3_rtw_ESSVR_alloc_hpw_pre  <- coxph(Surv(rtw_random, as.numeric(rtw_flg)) ~ alloc*hpw_pre, data=stroke)
summary(fit3_rtw_ESSVR_alloc_hpw_pre)

fit3_rtw_ESSVR_alloc_stroke_severity <- coxph(Surv(rtw_random, as.numeric(rtw_flg)) ~ alloc*stroke_severity, data=stroke)
summary(fit3_rtw_ESSVR_alloc_stroke_severity)

# output result of effect modification on ESSVR ~ Returning to work
fit3_eff_modi<- rbind(
  tidy(fit3_rtw_ESSVR_alloc_sex, conf.int = TRUE),
  tidy(fit3_rtw_ESSVR_alloc_age, conf.int = TRUE),
  tidy(fit3_rtw_ESSVR_alloc_work_status_pre, conf.int = TRUE),
  tidy(fit3_rtw_ESSVR_alloc_hpw_pre, conf.int = TRUE),
  tidy(fit3_rtw_ESSVR_alloc_stroke_severity, conf.int = TRUE)
  )
write.csv(fit3_eff_modi,"fit3_eff_modi.csv")

# AIM 2-2: Effect of ESSVR on returning to work adjusted for age, stroke severity and work_status_pre
fit3_rtw_ESSVR_alloc_adjusted <- coxph(Surv(rtw_random, as.numeric(rtw_flg)) ~ alloc+stroke_severity+work_status_pre+age, data=stroke)
summary(fit3_rtw_ESSVR_alloc_adjusted)


# AIM 2-2: Whether ESSVR affected return to work post-stroke time split

# define a function to save typing on the calculations
lincomb <- function(parnames, model, eform=TRUE) {
  npars <- length(parnames)
  cmat <- matrix(1, ncol=npars, nrow=npars)
  cmat[upper.tri(cmat)] <- 0
  b <- coef(model)[parnames]
  V <- vcov(model)[parnames, parnames]
  est <- as.vector(cmat %*% b)
  est <- est[npars]
  se <- sqrt(diag(cmat %*% V %*% t(cmat)))
  se <- se[npars]
  lincom <- c(exp(est), exp(est - 1.96*se), exp(est + 1.96*se), est/se, 2*pnorm(est/se, lower.tail = FALSE) )
  names(lincom) <- c("Estimate", "Lower.CI", "Upper.CI", "z value", "P value")
  if (eform) {
    lincom <- (lincom)
  }
  lincom 
}

fit3_rtw_ESSVR_timesplit <- coxph(Surv(tstart, rtw_random, as.numeric(rtw_flg)) ~ alloc * ftime_group, data=stroke_splite)
summary(fit3_rtw_ESSVR_timesplit)

lincomb(c("allocESSVR", "allocESSVR:ftime_group2"), fit3_rtw_ESSVR_timesplit)

# --------------------------------------------------------------------
# AIM 3: Whether the program affected health score (subgroup by sex) linear regression model
fit4_health_age <- lm(health_score ~ age, data=stroke)
fit4_health_sex <- lm(health_score ~ sex, data=stroke)
fit4_health_region <- lm(health_score ~ region, data=stroke)
fit4_health_work_status <- lm(health_score ~ work_status_pre, data=stroke)
fit4_health_hpw_pre <- lm(health_score ~ hpw_pre, data=stroke)
fit4_health_stroke_severity <- lm(health_score ~ stroke_severity, data=stroke)
fit4_health_alloc <- lm(health_score ~ alloc, data=stroke)

# summary of simple linear regression model
fit4_lm_summary <- rbind(tidy(fit4_health_sex, conf.int=TRUE)[2,],# significant
                         tidy(fit4_health_age, conf.int=TRUE)[2,],# significant
                         tidy(fit4_health_region, conf.int=TRUE)[2:4,],#significant
                         tidy(fit4_health_work_status, conf.int=TRUE)[2:5, ],#significant
                         tidy(fit4_health_hpw_pre, conf.int=TRUE)[2,], # not significant
                         tidy(fit4_health_stroke_severity, conf.int=TRUE)[2:3,], # significant
                         tidy(fit4_health_alloc, conf.int=TRUE)[2,] # not significant
                         )
write.csv(fit4_lm_summary, "lm_summary.csv")

# sub-group analysis on sex & measurement on effect modification
fit4_health_alloc_sex <- lm(health_score ~ alloc*sex, data=stroke)
summary(fit4_health_alloc_sex)

fit4_health_alloc_age <- lm(health_score ~ alloc*age, data=stroke)
summary(fit4_health_alloc_age)

fit4_health_alloc_region <- lm(health_score ~ alloc*region, data=stroke)
summary(fit4_health_alloc_region)

fit4_health_alloc_work_status_pre <- lm(health_score ~ alloc*work_status_pre, data=stroke)
summary(fit4_health_alloc_work_status_pre)

fit4_health_alloc_hpw_pre <- lm(health_score ~ alloc*hpw_pre, data=stroke)
summary(fit4_health_alloc_region)

fit4_health_alloc_stroke_severity <- lm(health_score ~ alloc*stroke_severity, data=stroke)
summary(fit4_health_alloc_stroke_severity)

fit4_modification_summary <- rbind(tidy(fit4_health_alloc_sex, conf.int=TRUE)[4,],# significant
                         tidy(fit4_health_alloc_age, conf.int=TRUE)[4,],# significant
                         tidy(fit4_health_alloc_region, conf.int=TRUE)[6:8,],#significant
                         tidy(fit4_health_alloc_work_status_pre, conf.int=TRUE)[7:10,],#significant
                         tidy(fit4_health_alloc_hpw_pre, conf.int=TRUE)[4,], # not significant
                         tidy(fit4_health_alloc_stroke_severity, conf.int=TRUE)[5:6,] # significant
                         )

write.csv(fit4_modification_summary, "lm_modification_summary.csv")

# multiple regression model adjusted by significant factors 
fit4_health_alloc_multi <- lm(health_score ~ alloc+sex+age+region+work_status_pre+stroke_severity, data=stroke)
summary(fit4_health_alloc_multi)
fit4_mlm_summary <- tidy(fit4_health_alloc_multi, conf.int=TRUE)
write.csv(fit4_mlm_summary, "fit4_mlm_summary.csv")
plot(fit4_health_alloc_multi, 1:2)

# sensitivity analysis ------------------------------------------------

stroke_sensi <- stroke %>% filter(essvr_complete_flg == 0 | is.na(essvr_complete_flg))
stroke_splite_sensi <- stroke_splite %>% filter(essvr_complete_flg == 0 | is.na(essvr_complete_flg))

fit3_rtw_ESSVR_sensi_timesplit <- coxph(Surv(tstart, rtw_random, as.numeric(rtw_flg)) ~ alloc * ftime_group, stroke_splite_sensi)
summary(fit3_rtw_ESSVR_sensi_timesplit)

lincomb(c("allocESSVR", "allocESSVR:ftime_group2"), fit3_rtw_ESSVR_sensi_timesplit)
