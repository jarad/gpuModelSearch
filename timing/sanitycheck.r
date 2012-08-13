## Sanity check on Nettleton's idea:

n = 10
rho = .5
x1 = rnorm(n)
x2 = rnorm(n,rho*x1)
y = rnorm(n,x1+2*x2) # truth: y = x1 + 2*x2 + e


y.x1 <- lm(y~x1-1)  
x2.x1 <- lm(x2~x1-1)

added.variable <- summary(lm(resid(y.x1)~resid(x2.x1)-1))

full <- summary(lm(y~x1+x2-1))

all.equal(coef(added.variable)[1],coef(full)[2,1])

all.equal(resid(added.variable), resid(full))


