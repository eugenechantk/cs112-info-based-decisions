library(foreign)
library(MASS)
library(Matching)
library(rgenoud)

set.seed(123)
"
######################################
Loading data set and editing the data
######################################
"
nsw_dw <- read_dta("nsw_dw.dta")
names(nsw_dw) <- c("data_id", "treat", "age", "educ", "black", "hisp", "married","nodgr","re74","re75","re78")

#transforming the dataset's data type
nsw_dw$data_id <- as.factor(nsw_dw$data_id)
nsw_dw$treat <- as.integer(nsw_dw$treat)
nsw_dw$age <- as.integer(nsw_dw$age)
nsw_dw$educ <- as.integer(nsw_dw$educ)
nsw_dw$black <- as.integer(nsw_dw$black)
nsw_dw$hisp <- as.integer(nsw_dw$hisp)
nsw_dw$married <- as.integer(nsw_dw$married)
nsw_dw$nodgr <- as.integer(nsw_dw$nodgr)
nsw_dw$re74 <- as.numeric(nsw_dw$re74)
nsw_dw$re75 <- as.numeric(nsw_dw$re75)
nsw_dw$re78 <- as.numeric(nsw_dw$re78)


"
##########################################################################
Running simple analysis on treatment effects on real income of 1978 (re78)
##########################################################################
"
#Point estimate of the difference in mean for nsw_dw dataset
pe.T = nsw_dw$treat
pe.Y = nsw_dw$re78
pe.mean.treated <- mean(pe.Y[pe.T==1])
pe.mean.untreated <- mean(pe.Y[pe.T==0])
pe.mean.treated  #mean for treated observation: 6349.144
pe.mean.untreated   #mean for untreated observation: 4554.801

#Simple, univariate linear regression to show treatment effect's statistical significance
lm.treat1 = lm(pe.Y~pe.T)
summary(lm.treat1)   #Treatment effect: 1794.3; treat coefficient has p-value of 0.00479, which is statistically significant
confint(lm.treat1, data=nsw_dw)  #95% CI of Tr effect: 550.5745, 3038.110


"
##########################################################
Replacing control units in nsw_dw with cps_control dataset
##########################################################
"
cps.control <- read_dta("cps_controls.dta")

#Formating the cps_control dataset to match nsw_dw dataset
names(cps.control) <- c("data_id", "treat", "age", "educ", "black", "hisp", "married","nodgr","re74","re75","re78")
cps.control$data_id <- as.factor(cps.control$data_id)
cps.control$treat <- as.integer(cps.control$treat)
cps.control$age <- as.integer(cps.control$age)
cps.control$educ <- as.integer(cps.control$educ)
cps.control$black <- as.integer(cps.control$black)
cps.control$hisp <- as.integer(cps.control$hisp)
cps.control$married <- as.integer(cps.control$married)
cps.control$nodgr <- as.integer(cps.control$nodgr)
cps.control$re74 <- as.numeric(cps.control$re74)
cps.control$re75 <- as.numeric(cps.control$re75)
cps.control$re78 <- as.numeric(cps.control$re78)

#Replacing cps_control to the control units in nsw_dw by creating a new dataset
new.nsw <- subset(nsw_dw,treat==1)
new.nsw <- rbind(new.nsw,cps.control)
attach(new.nsw)

#Simple point estimate and regression analysis on new dataset (185 treated, 15992 control)
T = treat
Y = re78
new.mean.treated <- mean(Y[T==1])
new.mean.untreated <- mean(Y[T==0])
new.mean.treated  #mean for treated observation: 6349.144
new.mean.untreated   #mean for untreated observation: 14846.66

#Simple, univariate linear regression to show treatment effect's statistical significance
lm.treat2 = lm(Y~T)
summary(lm.treat2)   #Treatment effect: -8497; treat coefficient has p-value of <0.0001, which is statistically significant
confint(lm.treat2)  #95% CI of Tr effect: -9893.156, -7101.877


"
#######################################
Matching 1.1: Propensity Score Matching
#######################################
"
#Match by propensity scores
pscore.glm <- glm(T~age+educ+black+hisp+married+nodgr+re74+re75,family = binomial(logit),data=new.nsw)  #construct propensity score
pscore = pscore.glm$fitted

pscore.match <- Match(Y=Y,Tr=T,X=pscore,M=1)  #Match based on the propensity score model above
summary(pscore.match)  #Est treatment effect: 1440.4 (p=0.0789); Due to large sample size, Tr. effect is statistically insignificant

match1.data <- pscore.match$mdata
with(match1.data, t.test(match1.data$Y~ match1.data$Tr))  #95% CI for Tr. effect: -1128.7872, -520.3607

"Check match balance"
match1.mb <- MatchBalance(T~pscore,data=new.nsw,match.out = pscore.match)
summary(match1.mb)

"Analysis on Match done by propensity score (in KS p-values):
The possible explanation for the statistically insignificant results is that the propensity score model is 
bias towards certain variables that are highly differentiated based on the p-value in the pscore model
 (e.g. black, hisp, married, nodgr and re75). Therefore, it fails to account for the difference in other covariates,
and lead to inexact matches."


"
###############################################################
Matching 1.2: Propensity Score Matching (include all covariate)
###############################################################
"
#Matching (including propensity score)
match.covariates <- cbind(age,educ,black,hisp,married,nodgr,re74,re75,pscore)
pscore.match2 <- Match(Y=Y,Tr=T,X=match.covariates,M=1)
summary(pscore.match2)  #Est. treatment effect: 1863.4 (p=0.0367); Tr. effect is statistically significant

match2.data <- pscore.match2$mdata
with(match2.data, t.test(match2.data$Y~ match2.data$Tr))  #95% CI for Tr. effect: -2191.0585, -666.3115

"Check for match balance"
match2.mb <- MatchBalance(T~age+educ+black+hisp+married+nodgr+re74+re75+pscore,data=new.nsw,match.out = pscore.match2)
summary(match2.mb)

"Analysis on match balance when matching with all covariates, including propensity score:
Covariates that matched well (in KS p-value): educ(p=0.998), black(p=1), hisp(p=1), married(p=1), nodgr(p=1), re75(p=0.912), pscore(p=0.982)
Covariates that didn't match well (in KS p-value): age (p=0.04), re74 (p<0.001)

Possible explanationn for poorly matched variables:
age and re74 both have a high p-value when constructing the propensity score model, indicating that the p-score will not
reflect the difference in age and re74 between treatment and control observations. Therefore, after matching, the difference
remains to be seen.
"


"
######################################################
Matching 2.1: Genetic Matching using Propensity Scores
######################################################
"
genmatch1.weights <- GenMatch(T,pscore,M = 1,pop.size = 10)
genmatch1 <- Match(Y,T,X=pscore,M=1,Weight.matrix = genmatch1.weights)
summary(genmatch1)  #Est. treatment effect: 1595.4 (p=0.05062); marginally significant given the large sample size

genmatch1.data <- genmatch1$mdata
with(genmatch1.data, t.test(genmatch1.data$Y~ genmatch1.data$Tr))  #95% CI for Tr. effect: -2556.222, -1060.703

"Check match balance"
genmatch1.mb <- MatchBalance(T~pscore,data=new.nsw,match.out = genmatch1)
summary(genmatch1.mb)

"Analysis on genetic matching on p-score:
Genetic matching delivered a better match compared to generic p-score matching. KS p-value for genetic matching is 0.15,
whereas generic matching is <0.001. However, the matching is not entirely exact.

The possible explanation for a better match compared to generic p-score matching is that the added weights to the p-score covariate 
increases the difference between some relatively dissimilar control observations that was originally tied. We can see that from the 
decreased number of matched observations (unnweighted). The number decreased from 4874 to 767. That means the treated observations 
are matched with less, but more similar control observations according to their p-scores
"


"
#########################################################################
Matching 2.2: Genetic Matching using all covariates and propensity scores
#########################################################################
"
genmatch2.covariates <- cbind(age,educ,black,hisp,married,nodgr,re74,re75,pscore)
genmatch2.weights <- GenMatch(T,genmatch2.covariates,M = 1,pop.size = 10,print.level = 0)
genmatch2 <- Match(Y,T,genmatch2.covariates,M=1,Weight.matrix = genmatch2.weights)
summary(genmatch2)  #Est. treatment effect: 1835.4 (p=0.052); marginally significant result

genmatch2.data <- genmatch2$mdata
with(genmatch2.data, t.test(genmatch2.data$Y~ genmatch2.data$Tr))  #95% CI for Tr. effects: -3189.598, -595.969

"Check match balance"
genmatch2.mb <- MatchBalance(T~age+educ+black+hisp+married+nodgr+re74+re75+pscore,data=new.nsw,match.out = genmatch2)
summary(genmatch2.mb)

"Analysis on genetic matching on all covariates and p-score:
Similar to the generic matching on all covariates and p-score, the covariates that matched well previously is also matched
well in the genetic matching. age's matching slightly improved (p=0.04 -> p=0.318), and so is re74 (p<0.001 -> p=0.036).
"


"
######################################
Extension: optimizing genetic matching
######################################
"
"Challenge: How to increase the balance for re74 (p=0.036) and age (p=0.318)"

exact_array = matrix(0,1,12)
exact_array[7] = 1
genmatch.opt.covariates = cbind(age,educ,black,hisp,married,nodgr,re74,re75,pscore,I(age*educ),I(age*married),I(age*nodgr))
genmatch.opt.weights <- GenMatch(T,genmatch.opt.covariates,M = 1,pop.size = 80,wait.generations = 10, caliper=0.25,print.level=0,exact = exact_array)
genmatch.opt <- Match(Y,T,genmatch.opt.covariates,M = 1,Weight.matrix = genmatch.opt.weights,caliper = 0.25,exact = exact_array)
summary(genmatch.opt)

"Check match balance"
genmatch.opt.mb <- MatchBalance(T~age+educ+black+hisp+married+nodgr+re74+re75+pscore,data=new.nsw,match.out = genmatch.opt)

"
Results (in KS p-value, or T-test p-value for dummy variables):
            pre-match     genmatch      gennmatch (optimized)
age         <0.0001       0.318         0.454                       
educ        <0.0001       0.954         1.000
black       <0.0001       1.000         1.000
hisp        <0.0001       1.000         1.000
married     <0.0001       1.000         1.000
nodgr       <0.0001       1.000         1.000
re74        <0.0001       0.036         1.000
re75        <0.0001       0.970         0.896
pscore      <0.0001       0.706         0.886

Parameteres optimized:
* Matching Covariates (genmatch.opt.covariates): I added a number of interaction terms to consider for matching. Matching on those
  interaction terms improved the balance on age (I got to p~0.8 without exact matching re74). I also tried to match on re74 interaction 
  terms, but they did not yield significant improvement on the balance.

* Exact matching on re74: After attempting to tweak different parameters on the GenMatch and adding interaction terms, I failed to improve
  on the balance of re74. Therefore, I decided to do an exact match on re74, leading to p=1.000.

* Increase population size on each generation: Since I did exact matching, it would reduce the balance on other covariates. Therefore, I
  increased the pop.size so that the GenMatch can dig deeper to find better weights to match for other covariates after exact-matching
  for re74.

* Increase threshold to stop populating new generations from 4 to 10: Similar rationale to increasing pop.size. To compensate on the
  imbalance of other covariates because of exact matching on re74, I want the GenMatch to dig deeper and find better weights for other
  covariates.

* Lowering the threshold distance of matching to 0.25 stdev: Lowering threshold distance would yield us closer matches across all the
  covariates. However, the cost of it is we discarded 112 treated observations. Instead of 185 matched observations, we are left with
  73, which may affect how conclusive our analysis is.
"