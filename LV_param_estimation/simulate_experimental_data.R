
library(deSolve)

# differentials along the time domain
chemLV <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
        dA <- -k1*A*X
        dX <- k1*A*X - k2*X*Y
        dY <- k2*X*Y - k3*Y
        dB <- k3*Y
        return(list(c(dA,dX,dY,dB)))
    })
}


pars  <- c(k1 = 0.3,
           k2 = 3,
           k3 = 0.1)

yini  <- c(A = 1,
           X = 1e-1,
           Y = 1e-1,
           B = 1e-6)

times <- seq(0, 500, by = 10) # every 5 second
out.ori <- ode(yini, times, chemLV, pars, method = "lsoda")
plot(out.ori, type='o')

# saving original (as it is easier to fit)
write.table(out.ori, file="lv_nontransformed.txt", sep='\t', quote=F, row.names = F)


# now let's do transforms that are typical for acquired biological data
#===============================================================================
out <- as.matrix(out.ori)
rownames(out) <- out[,1]
out <- out[,-1]

# 1. adding log-normal noise
# Note, the noise characteristics (st.dev.) are different between the species.
out.noise <- sapply(c(0.1,0.15,0.2,0.25), rnorm, n = nrow(out), mean=0)
out <- out * exp(out.noise)

# 2. censorship
# whatever is less is 1% of max value of a particular specie is set to be NA
# dynamic range of bioanalytical approaches is not infinite
out <- apply(out, 2, function(x){ x[x < max(x)*0.01] <- NA; x})

# 3. log2-transform and zero center
out <- log2(out)
out <- sweep(out, 2, colMeans(out, na.rm = TRUE), FUN='-')
#===============================================================================


# viz the rez
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
x <- as.data.frame(out) %>%
    rownames_to_column('time') %>%
    gather(specie,value,-time) %>%
    mutate(time = as.numeric(time))
    
ggplot(x, aes(x=time, y=value, color=specie)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ specie)


# saving
as.data.frame(out) %>%
    rownames_to_column('time') %>%
    write.table(file='lv.txt', sep='\t', quote=F, row.names=FALSE)

