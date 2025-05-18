# Eurovision Contest Televoting Monte Carlo simulation
# www.overfitting.net
# https://www.overfitting.net/


library(MCMCpack)  # Dirichlet distribution
library(tiff)


NCOUNTRIES=26  # 26 participating countries
N=1000000
gamma=2.2

for (distr in c('unif', 'unifdirichlet', 'nonunifdirichlet')) {
    set.seed(100)  # reproducible results
    if (distr=='unifdirichlet') {
        alpha=rep(1, NCOUNTRIES)
        votesnorm=t(rdirichlet(N, alpha))  # uniform Dirichlet distribution
    } else if (distr=='nonunifdirichlet') {
        alphamatrix=runif(NCOUNTRIES*N)  # define alpha weights with unif distribution
        dim(alphamatrix)=c(N, NCOUNTRIES)
        sampleslist=lapply(1:nrow(alphamatrix), function(i) {  # ver NOTA*
            rdirichlet(1, alphamatrix[i, ])  # non-uniform Dirichlet distribution
            })
        votesnorm=t(do.call(rbind, sampleslist))
    } else {
        votes=runif(NCOUNTRIES*N)  # uniform distribution
        dim(votes)=c(NCOUNTRIES, N)  # redim as NCOUNTRIES x N matrix
        votesnorm=sweep(votes, 2, colSums(votes), FUN="/")  # normalize each simulation
    }
    
    votesplot=votesnorm
    dim(votesplot)=c(dim(votesnorm)[1]*100, dim(votesnorm)[2]/100)
    writeTIFF(votesplot^(1/gamma), paste0("voting_", distr, ".tif"), bits.per.sample=16)
    
    minvotesrequired=apply(votesnorm, 2, max)*100  # get the min needed voting in %
    minimo=min(minvotesrequired)
    mediana=median(minvotesrequired)
    media=mean(minvotesrequired)
    intervalo95=quantile(minvotesrequired, 0.95)
    maximo=max(minvotesrequired)
    png(paste0("histogram_", distr, ".png"), width=1024, height=400)
        hist(minvotesrequired, breaks=800, xlim=c(0,100),
             main=paste0('Distr. of min % voting needed to win Eurovision televoting ',
                         "using a '", distr, "'", ' statistical distribution\n',
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
    dev.off()
    
    print(paste0("'", distr, "' simulation:"))
    print(paste0("  Televoto orgánico: ", 100/26, "%"))
    print(paste0("  Televoto promedio: ", media, "%"))
    print(paste0("  Televoto campaña: ", media-100/26, "%"))
}


# *NOTA: sobre el funcionamiento de la función con alpha de tipo matriz

# rdirichlet(n=nrow(alpha_matrix), alpha=alpha_matrix)
# con alpha_matrix matriz de 3 filas × 4 columnas, la función interpreta cada fila
# como un conjunto de parámetros Dirichlet distinto, y lo que devuelve es una
# matriz de tamaño n × k, en este caso 3 × 12 porque aplana los resultados
# en fila única concatenada.

# Esto es un comportamiento interno de rdirichlet() que no está documentado
# con claridad, pero sucede así:
    
# Si haces rdirichlet(3, alpha_vector_de_largo_4), obtienes 3 filas × 4 columnas
# como esperas.

# Si haces rdirichlet(3, alpha_matrix_3x4), devuelve 3 muestras personalizadas,
# pero como las aplana horizontalmente, obtienes una matriz de 3 × (3×4) = 3 × 12.
