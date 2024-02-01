library(pracma)

fi <- c(1/2,9/28,477/3136,543/21952,9433/19668992,4462/3813049,146689/1927561216,7155/1927561216,2809/1927561216)
m <- length(fi)-2
sum(fi)
mu <- sum(fi*(0:(m+1)))
fi_1 <- fi
Fi_1 <- cumsum(fi_1) # P(X<=x); x = 0,1,...,m,m+1
Fi_line_1 <- 1 - Fi_1 # P(X>x); x = 0,1,...,m,m+1
alpha_1 <- Fi_line_1[2:(m+1)]/fi_1[1] #alpha_k; k= 1,...,m
theroots <- polyroots(c(1,-alpha_1))
z <- theroots[[1]]
n <- theroots[[2]]
#
m <- length(alpha_1)
valphas <- rev(cumsum(rev(alpha_1)))
Malphas <- matrix(0,m,m)
for(i in 1:m){Malphas[i,i] <- 1}
for(j in 1:(m-1)){
  Malphas[(j+1):m,j] <- -alpha_1[1:(m-j)]
}
v_psi <- solve(Malphas,valphas)
l <- length(n)
Z <- matrix(0,m,m)
i <- 0
for(k in 1:l){
  for(j in 1:n[k]){
    i <- i + 1
    Z[,i] <- (1:m)^(j-1)*z[k]^(1:m)  
  }
}
bkj <- solve(Z,v_psi)
#
z1 <- Re(z[4])
b1 <- Re(bkj[7])
#
u <- c(1,2,4,6,8,10)
#
approx_psi_1 <- b1*z1^u
approx_psi_1
#
approx_psi_2 <-v_psi[1]*(v_psi[2]/v_psi[1])^(u-1)
approx_psi_2
