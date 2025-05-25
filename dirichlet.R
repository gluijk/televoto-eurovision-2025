# Dirichlet distribution
# www.overfitting.net
# https://www.overfitting.net/2025/05/interiorizando-posteriori-la.html


library(MCMCpack)  # Dirichlet distribution
library(tiff)


NVAR=3  # number of variables for the Dirichlet distribution
N=1000000  # number of Monte Carlo simulations
BREAKS=500  # histograms resolution
gamma=1  # 2.2

colores=c()  # we define up to 6 colours
colores=c(colores, rgb(1, 0, 0, alpha=0.3))  # red
colores=c(colores, rgb(0, 1, 0, alpha=0.3))  # green
colores=c(colores, rgb(0, 0, 1, alpha=0.3))  # blue
colores=c(colores, rgb(0, 1, 1, alpha=0.3))  # red inv
colores=c(colores, rgb(1, 0, 1, alpha=0.3))  # green inv
colores=c(colores, rgb(1, 1, 0, alpha=0.3))  # blue inv
coloressolid=sub("(.{6}).{2}$", "\\1FF", colores)


distr=c('unif', 'unifdirichlet_1_1_1',
        'nonunifdirichlet_0.9_0.9_0.9', 'nonunifdirichlet_1.5_1.5_1.5',
        'nonunifdirichlet_5_5_5', 'nonunifdirichlet_1_2_2',
        'nonunifdirichlet_2_4_8', 'nonunifdirichlet_random')

for (j in 1:length(distr)) {
    set.seed(100)  # reproducible results
    if (distr[j]=='unif') {
        votes=runif(N*NVAR)  # uniform
        dim(votes)=c(N, NVAR)  # redim as N x NVAR matrix
        votesnorm=sweep(votes, 1, rowSums(votes), FUN="/")  # normalize sum=1
    } else if (distr[j]=='unifdirichlet_1_1_1') {
        alpha=rep(1, NVAR)
        votesnorm=rdirichlet(N, alpha)  # non-uniform symmetric Dirichlet
    } else if (distr[j]=='nonunifdirichlet_0.9_0.9_0.9') {
        alpha=rep(0.5, NVAR)
        votesnorm=rdirichlet(N, alpha)  # non-uniform symmetric Dirichlet
    } else if (distr[j]=='nonunifdirichlet_1.5_1.5_1.5') {
        alpha=rep(1.5, NVAR)
        votesnorm=rdirichlet(N, alpha)  # non-uniform symmetric Dirichlet
    } else if (distr[j]=='nonunifdirichlet_5_5_5') {
        alpha=rep(5, NVAR)
        votesnorm=rdirichlet(N, alpha)  # non-uniform asymmetric Dirichlet
    } else if (distr[j]=='nonunifdirichlet_1_2_2') {
        alpha=c(1, 2, 2)
        votesnorm=rdirichlet(N, alpha)  # non-uniform asymmetric Dirichlet
    } else if (distr[j]=='nonunifdirichlet_2_4_8') {
        alpha=c(2, 4, 8)
        votesnorm=rdirichlet(N, alpha)  # non-uniform Dirichlet
    } else if (distr[j]=='nonunifdirichlet_random') {
        alphamatrix=runif(N*NVAR)  # define alpha weights with unif
        dim(alphamatrix)=c(N, NVAR)
        sampleslist=lapply(1:nrow(alphamatrix), function(i) {  # ver NOTA*
            rdirichlet(1, alphamatrix[i, ])  # non-uniform random Dirichlet
            })
        votesnorm=do.call(rbind, sampleslist)
    }

    mediana=c()
    medianas="partial medians="
    MAXIMO=0
    for (i in 1:NVAR) {
        mediana=c(mediana, median(votesnorm[,i]))
        medianas=paste0(medianas, ifelse(i>1, ' / ',''), round(mediana[i], 2))
        MAXTMP=max(hist(votesnorm[,i], breaks=BREAKS, plot=FALSE)$counts)
        if (MAXTMP>MAXIMO) MAXIMO=MAXTMP
    }
    media=mean(votesnorm)  # 1/3 by definition
    
    name=paste0(j, "_", distr[j])

    # Plot histograms of min % voting needed
    png(paste0("histogram", name, ".png"), width=800, height=400)
        hist(votesnorm[,1], breaks=BREAKS, col=colores[1],
             xlim=c(0,1), ylim=c(0,MAXIMO),
             main=paste0("'", distr[j], "'", ' distribution\n',
                         '(',
                         medianas, ', ',
                         'global mean=', round(media,2),
                         ')'
                         ),
            xlab='x, y, z values', border=NA)
        abline(v=mediana[1], col=coloressolid[1], lty="dotted")

        for (i in 2:NVAR) {
            hist(votesnorm[,i], breaks=BREAKS, col=colores[i],
                               add=TRUE, border=NA)
            abline(v=mediana[i], col=coloressolid[i], lty="dotted")
        }
        abline(v=media)
    dev.off()
    
    # Draw Dirichlet Triangle
    if (NVAR==3) {  # only for 3 variables the triangle makes sense
        df=as.data.frame(votesnorm)
        colnames(df)=c('x', 'y', 'z')
        df$u = (-df$x + df$y +          1) / 2  # 2D projection of the x+y+z=1 plane
        df$v = (-df$x - df$y + 2*df$z + 1) / (2*3^0.5)
        png(paste0("triangle", name, ".png"), width=800, height=800)
        plot(df$u, df$v,
             main=paste0("'", distr[j], "'", ' statistical distribution'),
             xlim=c(0,1), ylim=c(0,0.9),
             pch=19, cex=0.25, asp=1,
             col=rgb(0, 0, 1, alpha=0.01),
             xlab="", ylab="")
        lines(c(0, 1, 1/2, 0), c(0, 0, 3^0.5/2, 0), lty="dotted")
        dev.off()
    }
    
    # Render graphically the voting distribution
    votesplot=votesnorm
    dim(votesplot)=c(dim(votesnorm)[1]/1000, dim(votesnorm)[2]*1000)
    writeTIFF(votesplot^(1/gamma), paste0("voting", name, ".tif"),
              bits.per.sample=16)
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
