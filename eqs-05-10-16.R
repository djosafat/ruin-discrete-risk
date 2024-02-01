library(pracma)
###############################
p_prime <- function(x,alpha){
  m <- length(alpha)
  aux <- 0
  for(k in 1:m){
    aux <- aux + alpha[k]*(m-k)*x^(m-k-1)
  }
  return(m*x^(m-1)-aux)
}
###############################
eq_5 <- function(u,fi,mu){
        Fi <- cumsum(fi) # P(X<=x); x = 0,1,...,m,m+1
        Fi_line <- 1 - Fi # P(X>x); x = 0,1,...,m,m+1
        m <- length(fi) - 2  
        if(u == 0){
                return(mu)
        }else if(u == 1){
                return(sum(Fi_line[2:(m+2)])/fi[1])
        }else{
                aux1 <- 0
                if(u<=(m+1)){
                        for(k in 1:(u-1)){
                                aux1 <- aux1 + Fi_line[k+1]*eq_5(u-k,fi,mu)
                        }
                        aux2 <- sum(Fi_line[(u+1):(m+2)])
                }else{
                        for(k in 1:(m+1)){
                                aux1 <- aux1 + Fi_line[k+1]*eq_5(u-k,fi,mu)
                        }
                        aux2 <- 0
                        
                }
                return((aux1+aux2)/fi[1])
        }
}
###############################
eq_10 <- function(u,mu,fi,alpha,z){
    m <- length(alpha)  
    aux <- (1-mu)/fi[1]
    sk <- (z[1]^(m+u-1))/((1-z[1])*(p_prime(z[1],alpha)))
    if(length(z)>1){
      for(k in 2:m){
        num <- z[k]^(m+u-1)
        den <- (1-z[k])*(p_prime(z[k],alpha))
        sk <- sk + num/den
      }
    }
    return(aux*sk)
}
################################
eq_16 <- function(u,mu,fi,alpha,z,n){
  m <- length(alpha)
  valphas <- rev(cumsum(rev(alpha)))
  Malphas <- matrix(0,m,m)
  for(i in 1:m){Malphas[i,i] <- 1}
  for(j in 1:(m-1)){
    Malphas[(j+1):m,j] <- -alpha[1:(m-j)]
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
  aux <- 0
  i <- 0
  for(k in 1:l){
    for(j in 1:n[k]){
      i <- i + 1
      aux <- aux + bkj[i]*u^(j-1)*z[k]^u
    }
  }
  return(aux)
}
################################
psi <- function(u,fi){
    m <- length(fi)-2
    mu <- sum(0:(m+1)*(fi))
    if(u == 0){
        return(mu)
    }else{
        Fi <- cumsum(fi) # P(X<=x); x = 0,1,...,m,m+1
        Fi_line <- 1 - Fi # P(X>x); x = 0,1,...,m,m+1
        alpha <- Fi_line[2:(m+1)]/fi[1] #alpha_k; k= 1,...,m
        theroots <- polyroots(c(1,-alpha))
        z <- theroots[[1]]
        n <- theroots[[2]]
        if(length(z)==m){  #simple roots
            eq_10(u,mu,fi,alpha,z)    
        }else{  # multiplicities > 1
            eq_16(u,mu,fi,alpha,z,n)
        }
    }
}

#####################
#Example section 6.1
#####################
fi <- c(7/8,0,0,0,0,0,0,1/8)
m <- length(fi)-2
sum(fi)
mu <- sum(fi*(0:(m+1))) # 7/8
fi_1 <- fi
Fi_1 <- cumsum(fi_1) # P(X<=x); x = 0,1,...,m,m+1
Fi_line_1 <- 1 - Fi_1 # P(X>x); x = 0,1,...,m,m+1
alpha_1 <- Fi_line_1[2:(m+1)]/fi_1[1] #alpha_k; k= 1,...,m
theroots_1 <- polyroots(c(1,-alpha_1))
z_1 <- theroots_1[[1]] #roots
n_1 <- theroots_1[[2]] #multiplicities
############
v <- c(1,12,24,36,48,60)
N <- length(v)
psi_v <- c()
for(u in v){
  psi_v <- c(psi_v,psi(u,fi))
}
Re(psi_v) # exact ruin probability using eq (10) or (16)
plot(v,Re(psi_v),type = 'l')
#############################################
psi_classic_v <- c()
for(u in v){
  psi_classic_v <- c(psi_classic_v,eq_5(u,fi,mu))
}
psi_classic_v # exact ruin probability using eq (5) very slow
#############################################

#####################
#Example section 6.2
#####################
fi <- c(1/2,9/28,477/3136,543/21952,9433/19668992,4462/3813049,146689/1927561216,7155/1927561216,2809/1927561216)
m <- length(fi)-2
sum(fi)
mu <- sum(fi*(0:(m+1)))
fi_1 <- fi
Fi_1 <- cumsum(fi_1) # P(X<=x); x = 0,1,...,m,m+1
Fi_line_1 <- 1 - Fi_1 # P(X>x); x = 0,1,...,m,m+1
alpha_1 <- Fi_line_1[2:(m+1)]/fi_1[1] #alpha_k; k= 1,...,m
theroots_1 <- polyroots(c(1,-alpha_1))
z_1 <- theroots_1[[1]] #roots
n_1 <- theroots_1[[2]] #multiplicities
############
uu <- 10
############
v <- c(1,2,4,6,8,10)
N <- length(v)
psi_v <- c()
for(u in v){
        psi_v <- c(psi_v,psi(u,fi))
}
Re(psi_v) # exact ruin probability using eq (10) or (16)
plot(v,Re(psi_v),type = 'l') 
#############################################
psi_classic_v <- c()
for(u in v){
        psi_classic_v <- c(psi_classic_v,eq_5(u,fi,mu))
}
psi_classic_v # exact ruin probability using eq (5)
#############################################
fi_2 <- fi
epsilon <- 10^(-3)
fi_2[1] <- fi_2[1]+epsilon
fi_2[2] <- fi_2[2]-epsilon
sum(fi_2)
Fi_2 <- cumsum(fi_2) # P(X<=x); x = 0,1,...,m,m+1
Fi_line_2 <- 1 - Fi_2 # P(X>x); x = 0,1,...,m,m+1
alpha_2 <- Fi_line_2[2:(m+1)]/fi_2[1] #alpha_k; k= 1,...,m
sum(alpha_1)
sum(alpha_2)
theroots_2 <- polyroots(c(1,-alpha_2))
theroots_1
theroots_2
abs(theroots_1)
abs(theroots_2)
############
psi_2_v <- c()
for(u in v){
        psi_2_v <- c(psi_2_v,psi(u,fi_2))
}
Re(psi_2_v) # approximation epsilon simple roots
abserror <- abs(Re(psi_2_v)-Re(psi_v))/Re(psi_v)
abserror
plot(v,abserror,col='blue',type='l')
