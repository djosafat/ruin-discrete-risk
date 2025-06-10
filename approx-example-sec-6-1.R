library(pracma)
##########################################
fi <- c(7/8,0,0,0,0,0,0,1/8)
m <- length(fi)-2
sum(fi)
mu <- sum(fi*(0:(m+1))) # 7/8
fi_1 <- fi
Fi_1 <- cumsum(fi_1) # P(X<=x); x = 0,1,...,m,m+1
Fi_line_1 <- 1 - Fi_1 # P(X>x); x = 0,1,...,m,m+1
alpha_1 <- Fi_line_1[2:(m+1)]/fi_1[1] #alpha_k; k= 1,...,m
theroots_1 <- polyroots(c(1,-alpha_1))
z <- theroots_1[[1]] #roots
l <- m
m <- length(alpha_1)
valphas <- rev(cumsum(rev(alpha_1)))
Malphas <- matrix(0,m,m)
for(i in 1:m){Malphas[i,i] <- 1}
for(j in 1:(m-1)){
  Malphas[(j+1):m,j] <- -alpha_1[1:(m-j)]
}
v_psi <- solve(Malphas,valphas)
Z <- matrix(0,m,m)
i <- 0
for(k in 1:l){
    i <- i + 1
    Z[,i] <- z[k]^(1:m)  
}
bkj <- solve(Z,v_psi)
#
z1 <- Re(z[6])
b1 <- Re(bkj[6])
#
u <- c(1,12,24,36,48,60)
#
approx_psi_1 <- b1*z1^u
approx_psi_1
#
psi1 <- v_psi[1]
psi2 <- v_psi[2]
approx_psi_2 <-psi1*(psi2/psi1)^(u-1)
approx_psi_2
