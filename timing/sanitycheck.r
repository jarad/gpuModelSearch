## Sanity check on Nettleton's idea:
## Full model is y = a + xb + zg + d
## Red model is  y = a + xb + e
## Attempt to capture the residuals of the full model by regressing the
## residuals of the reduced models on the new covariate, e.g. with
## ehat = zg + d or possibly ehat = a + zg + d

set.seed(154)
dat <- data.frame(y=rnorm(5, mean=5, sd=2), x=rnorm(5, mean=5, sd=2), z=rnorm(5, mean=5, sd=2))

xslr <- lm(y~x, data=dat)  ## simple linear regression: y = a + xb + e
xz <- lm(y~x+z, data=dat)  ## full model
dat$rx <- resid(xslr)

xz2noint <- lm(rx~z-1, data=dat) ## ehat = zg + d
xz2int <- lm(rx~z, data=dat)     ## ehat = a + zg + d

resids <- data.frame(xz=resid(xz), xz2noint=resid(xz2noint), xz2int=resid(xz2int))

## None of the resids match
resids

SSR <- rep(0,3)
for(i in 1:3)
    SSR[i] <- sum(resids[,i]^2)

## The sums of squared residuals also don't match
SSR
