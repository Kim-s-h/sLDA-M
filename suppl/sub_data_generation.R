# ------------------ Data generation ---------------------- #
nDocs = 1000   # number of documents
K = 4         # number of topics
V = 400       # vocabulary size
nWords = 100  # document length
nG = 2        # number of groups

par.b = 5
var.cond = 2  # error variance


b1.coef = c(4, 8, 8, 4)  # regression coefficient for Group 1
b2.coef = c(3, 5, 5, 3)  # regression coefficient for Group 2


#beta
betas = matrix(1,ncol=V,nrow=K)
for (v in 1:100){
  betas[1,v] = par.b
}
for (v in 101:200){
  betas[2,v] = par.b
}
for (v in 201:300){
  betas[3,v] = par.b
}
for (v in 301:400){
  betas[4,v] = par.b
}


## Generate Topics
gam = matrix(NA,nrow=K, ncol=V)
for (k in 1:K){
  gam[k,] = rdirichlet(1,alpha=betas[k,])
}

# plot(gam[1,], type='l')
# points(gam[2,], type='l', col='tomato')
# points(gam[3,], type='l', col='steelblue')
# points(gam[4,], type='l', col='green4')


## Generate membership
mem = sample(c(1:nG),nDocs, replace = T)

## Generate theta
theta = rdirichlet(nDocs, rep(1,K))

## Generate document length
N.D = rpois(nDocs,lambda=nWords) 

## Generate words
y_kv = matrix(0,nrow = K, ncol = V )
x_dk = matrix(0,nrow = nDocs, ncol = K)
Z_t = rep(NA, sum( N.D ) )
W = matrix(NA, nrow = sum(N.D), ncol = 2)
index = 0
for(d in 1:nDocs){
  for(n in 1:N.D[d] ){
    index = index + 1
    prob_Z = theta[d,]
    Z_t[ index ] = sample(1:K, 1, prob = prob_Z)
    prob_V = gam[Z_t[index],]
    W[index,1] = sample(1:V, 1, prob = prob_V )
    W[index,2] = d
    y_kv[ Z_t[index], W[index,1] ] = 1 + y_kv[ Z_t[index], W[index,1] ]
    x_dk[ W[index,2], Z_t[index] ] = 1 + x_dk[ W[index,2], Z_t[index] ]
  }
}

## Generate response
zbar = t(apply(x_dk, 1, function(x) x/sum(x)))

resp = c()
for(i in 1:nDocs){
  if(mem[i] == 1){
    resp[i] = rnorm(1, mean=t(b1.coef)%*%zbar[i,], sd=sqrt(var.cond))
  } else if(mem[i] == 2){
    resp[i] = rnorm(1, mean=t(b2.coef)%*%zbar[i,], sd=sqrt(var.cond))
  }
}


t.sig = var.cond
t.gam = gam
t.theta = theta
t.Z = Z_t
t.coef1 = b1.coef
t.coef2 = b2.coef

