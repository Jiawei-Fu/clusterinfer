# clusterinfer
R package for statistical inference involving small cluster sizes

The package is used to improve the statistical inference involving  small cluster sizes. Instead of calculating robust cluster variance-covariance  matirx or bootstrap-se, the functions calculating p values based on  wild bootstrap-t for OLS models and score based wild bootstrap-t for MLE models. The refined p value significantly reduces the probability of over-rejection of  the null hypothesis. 

```{r}
install.packages("devtools")

```
