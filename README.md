# clusterinfer
R package for statistical inference involving small cluster sizes (smaller than 50)

Author: Jiawei Fu

Contributor: Daniel Stegmueller

The package is used to improve the statistical inference involving  small cluster sizes. Instead of calculating robust cluster variance-covariance  matirx or bootstrap-se, the functions calculating p values based on  wild bootstrap-t for OLS models and score based wild bootstrap-t for MLE models. The refined p value significantly reduces the probability of over-rejection of  the null hypothesis. 

To install and use the latest version of the package, try the following codes:
```r
install.packages("devtools")
devtools::install_github("Jiawei-Fu/clusterinfer")
library(clusterinfer)
```
Some of users will receive warnings about time zone, it will not affect the installation.

There are two functions and one data set in the package now.

The function wildboott is designed to OLS linear models.

The function is used to get accurate and precise inference from small cluster sizes. Through  (restricted) wild clustered bootstrap-t to get p-values for OLS models, the function considerably reduces the probability of over rejection of the null hypotheses.The default distribution used to reconstruct the residuals is Rademacher. When cluster size is smaller than 12, the algorithm will suggest using "six - point" distribution instead.

Example:
```r
data(wv6_equ)  # input data
wv6_equ <- as.data.frame(wv6_equ)
model_equ <- lm(income_equ ~ income + age + gender, data = wv6_equ) # linear model
summary(model_equ)
wildboott(model_equ, ~country, R = 250, seed = 320) #find original p values are over estimated
```
The function swildboott is designed to MLE nonlinear models.

The function is used to get accurate and precise inference from small cluster sizes. Through  (restricted) wild clustered bootstrap-t to get p-values for MLE models, the function considerably reduces the probability of over rejection of the null hypotheses. The distribution used to reconstruct the residuals is Mammen or Rademacher.

Examples:
```r
##### glm model #####
data(wv6_equ)  # input data
# create binomial dependent variable
wv6_bi <- as.data.frame(wv6_equ)
wv6_bi$income_equ0[wv6_bi$income_equ < 6] <- 0
wv6_bi$income_equ0[wv6_bi$income_equ > 5] <- 1
glm_equ <- glm(income_equ0 ~ income + age + gender, data = wv6_equ) # glm model
swildboott(glm_equ, ~country, R = 500) # find original p values are over estimated

##### polr model #####
data(wv6_equ)  # input data
require(MASS)
wv6_equ <- as.data.frame(wv6_equ)
polr_equ <- polr(factor(income_equ) ~ income + age + gender, data = wv6_equ) # polr model
swildboott(polr_equ, ~country, R = 500, type = "mammen") # find original p value is over estimated
```

Next version will be avaliable to more regreesion models and to deal with multi-way clustering and hybrids estimators.
If you have any comments, suggestions, or findings, appreciate you sending me the email to: jiawei.fu@duke.edu

