library(Matching)
library(boot)
library(tree)
library(randomForest)
data("lalonde")

"
Two multiple regression models:
nodegr.lm.fit: no-degree subgroup; ind-var:re78; dep-var:treat, age, black, hisp
degr.lm.fit: degree subgroup; ind-var:re78; dep-var:treat, age, black, hisp
"
nodegr.lm.fit=lm(re78~treat+age+black+hisp, data=subset(lalonde, nodegr==1))
degr.lm.fit=lm(re78~treat+age+black+hisp, data=subset(lalonde, nodegr==0))

"Finding the p-values of each dep-var in the regressions"
summary(nodegr.lm.fit)
summary(degr.lm.fit)

"Finding the confidence interval of each of the dep-var in the regressions"
confint(nodegr.lm.fit)
confint(degr.lm.fit)

###

"
Simple regression models:
inter.lm.fit: ind-var:re78; dep-var:treat*nodegr
"
inter.lm.fit=lm(re78~treat*nodegr, data=lalonde)

"Coefficient and p-values of the interaction simple regression"
summary(inter.lm.fit)

"
BOOTSTRAPPING A SIMPLE REGRESSION MODEL TO OBTAIN STD ERR
"
"Create a bootstrap function, with the simple regression model as the output"
boot.fn = function(data,index){
  return(coef(lm(re78~treat*nodegr, data=data,subset = index)))
}
boot.fn(lalonde,1:445)
set.seed(1)
"Bootstrap 1000 times"
boot.inter <- boot(lalonde,boot.fn,2000)
#boot.inter
"Calculate BASIC confidence interval with bootstrap"
boot.inter.ci <- boot.ci(boot.inter,conf = 0.95)
#boot.inter.ci
ci.t = boot.inter.ci$basic[,c(4,5)]

hist(boot.inter$t[,2], main = 'Confidence Interval for Treatment Effect', xlab = 'RSquared',
     col = 'grey')
abline(v = ci.t, col = 'red',lwd=2)
abline(v=c(-172,2564),col=c('blue','blue'))
abline(v=c(105,6194),col=c('green','green'))
abline(v=mean(boot.inter$t[,2]),col='red',lty=2,lwd=2)
legend("topright", inset = 0.01, legend=c("Non-diploma", "Diploma", "Treat*Nodegr", "Est. Effect"),
       col=c("blue", "green","red","red"),lty=c(1,1,1,2),lwd=c(1,1,2,2), cex=0.8)

###

"Create u78 (unemployment in 1978) in lalonde data set"
lalonde$u78[lalonde$re78<=0] <- 1
lalonde$u78[lalonde$re78>0] <- 0
lalonde <- lalonde[c(1,2,3,4,5,6,7,8,9,10,11,13,12)]

"
Two multiple logistic regression models:
nodegr.lm.fit: no-degree subgroup; ind-var:u78; dep-var:treat, age, black, hisp
degr.lm.fit: degree subgroup; ind-var:u78; dep-var:treat, age, black, hisp
"
nodegr.glm.fit=glm(u78~treat+age+black+hisp,family=binomial,data=subset(lalonde,nodegr==1))
degr.glm.fit=glm(u78~treat+age+black+hisp,family=binomial,data=subset(lalonde,nodegr==0))

"Finding the treatment coefficient in the logistic regressions"
summary(nodegr.glm.fit)$coef
summary(degr.glm.fit)$coef

"Calculating the odds of the treatment --> treatment effect"
exp(coef(nodegr.glm.fit))
exp(coef(degr.glm.fit))

"Calculating confidence interval using standard error of coefficient in the logistic regressions"
exp(confint.default(nodegr.glm.fit))
exp(confint.default(degr.glm.fit))

###

"Random Forest for non-diploma holders; 3 predictors considered; 15 trees to minimize error"
ranfor.u78.nodegr=randomForest(u78~.-re78,data=subset(lalonde,nodegr==1),mtry=3,ntree=15)
#plot(ranfor.u78.nodegr) --> used to test for error
"Variable importance plot for random forest of non-diploma holders"
varImpPlot(ranfor.u78.nodegr, main="Variable Importance for Non-Diploma Holders")

"Random Forest for non-diploma holders; 3 predictors considered; 22 trees to minimize error"
ranfor.u78.degr=randomForest(u78~.-re78,data=subset(lalonde,nodegr==0),mtry=3,ntree=17)
#plot(ranfor.u78.degr) --> used to test for error
"Variable importance plot for random forest of diploma holders "
varImpPlot(ranfor.u78.degr, main="Variable Importance for Diploma Holders")
