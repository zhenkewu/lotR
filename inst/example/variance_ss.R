prob = 0.2
mu1 = 2
mu0 = 0

sigma1 = 3
sigma0 = 1

N = 1000000
s = rbinom(N,1,prob)
y <- rep(NA,length(s))
for (i in seq_along(s)){
      y[i] <- rnorm(1,mu1,sqrt(sigma1))*s[i] +
        rnorm(1,mu0,sqrt(sigma0))*(1-s[i])
}

var(s*y)
prob*(sigma1+(1-prob)*mu1^2)


mean(y^2)
prob*(sigma1+mu1^2) + (1-prob)*sigma0

mean(s*y^2)
prob*(sigma1+mu1^2)

mean(s*y)
prob*mu1
