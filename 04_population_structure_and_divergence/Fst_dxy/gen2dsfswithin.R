library(genio)
library(parallel)
getFst<-function(est){ # from https://raw.githubusercontent.com/ANGSD/angsd/master/R/fstFrom2dRHELLER.R
    N1<-nrow(est)-1
    N2<-ncol(est)-1
    # cat("N1: ",N1 ," N2: ",N2,"\n")
    est0<-est
    est0[1,1]<-0
    est0[N1+1,N2+1]<-0
    est0<-est0/sum(est0)

    aMat<<-matrix(NA,nrow=N1+1,ncol=N2+1)
    aMat.ss<<-matrix(NA,nrow=N1+1,ncol=N2+1)
    baMat<<-matrix(NA,nrow=N1+1,ncol=N2+1)
    for(a1 in 0:(N1))
        for(a2 in 0:(N2)){
            p1 <- a1/N1
            p2 <- a2/N2
            q1 <- 1 - p1
            q2 <- 1 - p2
            N <- (p1-p2)^2
            D <- p1*(1-p2)+p2*(1-p1)
            aMat[a1+1,a2+1]<<-N
            baMat[a1+1,a2+1]<<-D
            #sample size correction
            N.ss <- (p1-p2)^2-((p1*(1-p1))/(N1-1))-((p2*(1-p2))/(N2-1))
            aMat.ss[a1+1,a2+1]<<-N.ss
        }
    ## sample size corrected moment estimator
    ss <- sum(est0*aMat.ss,na.rm=T)/sum(est0*baMat,na.rm=T) 
    c(fstSS=ss)
}

plink <- read_plink("/projects/seqafrica/people/mdn487/Topi/goatBased/Goat_topi_variable_biallelic_noindels_polymorphic_goodSites")
fam <- as.data.frame(plink$fam)
popinfo <- fam$fam
indnames <- fam$id
colnames(plink$X) <- indnames
print("plink file loaded")

totsitesn <- 1065447319 # magic number used to scale dxy from good sites filter

X <- as.matrix(plink$X)
combs <- t(combn(indnames, 2))

getSFS <- function(index, geno) {
    sfs <- table(X[, combs[index, 1]], X[, combs[index, 2]])
    sfs[1, 1] <- totsitesn - (sum(sfs) - sfs[1, 1])
    fst <- getFst(sfs)
    dxy <- (2 * (sfs[2, 2] + sfs[2, 3] + sfs[3, 2] + sfs[1, 2] + sfs[2, 1]) + 4 * (sfs[1, 3] + sfs[3, 1])) / (4 * sum(sfs))
    cat(index, "\n")
    c(as.vector(sfs), unname(fst), dxy)
}
sfslist <- mclapply(1:nrow(combs), getSFS, geno = X, mc.preschedule = TRUE, mc.cores = 64)

mat <- do.call(rbind, sfslist)
mat <- cbind(combs, mat)
colnames(mat) <- c("A", "B", "A0B0", "A1B0", "A2B0", "A0B1", "A1B1", "A2B1", "A0B2", "A1B2", "A2B2", "fst", "dxy")
write.table(mat, file = paste0("./goat.topi.allsites.2dsfs.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

