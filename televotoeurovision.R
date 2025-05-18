# Eurovision Contest Televoting Monte Carlo simulation
# www.overfitting.net
# https://www.overfitting.net/


NCOUNTRIES=26  # 26 participating countries
N=10000000  # number of simulations
set.seed(100)  # reproducible results
votes=runif(NCOUNTRIES*N)
dim(votes)=c(NCOUNTRIES, N)  # redim as NCOUNTRIES x N matrix
votesnorm=sweep(votes, 2, colSums(votes), FUN="/")  # normalize each simulation
minvotesrequired=apply(votesnorm, 2, max)*100  # get the min needed voting in %

minimo=min(minvotesrequired)
mediana=median(minvotesrequired)
media=mean(minvotesrequired)
intervalo95=quantile(minvotesrequired, 0.95)
maximo=max(minvotesrequired)
hist(minvotesrequired, breaks=800, xlim=c(0,ceiling(maximo/10)*10),
     main=paste0('Distribution of min % of voting needed to win Eurovision televoting\n',
                 '(',
                 'min=', round(minimo,1), '%, ',
                 'median=', round(mediana,1), '%, ',
                 'mean=', round(media,1), '%, ',
                 'conflevel95=', round(intervalo95,1), '%, ',
                 'max=', round(maximo,1), '%',
                 ')'
                 ),
    xlab='% of voting')
abline(v=c(media, mediana), col='red')
abline(v=intervalo95, col='red', lty="dotted")
abline(v=c(minimo, maximo), col='gray', lty="dotted")

print(paste0("Televoto orgánico: ", 100/26, "%"))
print(paste0("Televoto promedio: ", media, "%"))
print(paste0("Televoto campaña: ", media-100/26, "%"))