library(quantmod)
library(moments) #skew and kurtosis
library(MASS) #fitdistr
library(metRology) #rt.scaled
library(rugarch)  #ugarspec , ugarfit


wilsh <- getSymbols("WILL5000IND",src="FRED",auto.assign = FALSE)
brit<-getSymbols("DEXUSUK",src="FRED",auto.assign = FALSE)
jap<-getSymbols("DEXJPUS",src="FRED",auto.assign = FALSE)
swiss <- getSymbols("DEXSZUS",src="FRED",auto.assign = FALSE)
aus <- getSymbols("DEXUSAL",src="FRED",auto.assign = FALSE)
gold <- getSymbols("GOLDPMGBD228NLBM",src="FRED",auto.assign = FALSE)


wilsh <-na.omit(wilsh)
wilsh <-wilsh["1979-12-31/2017-12-31"]

#gold <-gold["1979-12-31/2008-09-15"]
#gold <-gold["1979-12-31/1987-10-19"]

jap$DEXJPUS = 1/jap$DEXJPUS
swiss$DEXSZUS = 1/swiss$DEXSZUS
#names(wilsh) <-"TR"
logret <-diff(log(wilsh))[-1]     #daily log return
#other compounded returns
logret.w <-apply.weekly(logret,sum)
logret.m <-apply.monthly(logret,sum)
logret.q <- apply.quarterly(logret,sum)
logret.y <- apply.yearly(logret,sum)

#discrete return
ret <- exp(logret)-1

#Direct ES
es <- mu-sig*dnorm(qnorm(0.05,0,1),0,1)/0.05

#vector
rvec <-as.vector(logret)

#general terms
mu = mean(logret)
sig = sd(logret)

#for normal assumed Calculations
rvec <-rnorm(100000,mu,sig)
#without normal assumed
rvec <- sample(as.vector(logret),100000,replace=TRUE)

VaR <-quantile(rvec,alpha)
ES <-mean(rvec[rvec<VaR])
round(VaR,6)
round(ES,6)





#normality variables
round(skewness(rvec),2)
round(kurtosis(rvec),2)
jarque.test(rvec)



t.fit<-fitdistr(rvec,'t')
round(t.fit$estimate,6)

#direct student-t distributed for 1 day
alpha <-0.05
set.seed(123789)
rvec<-rt.scaled(100000,mean=t.fit$estimate[1],sd=t.fit$estimate[2],df=t.fit$estimate[3])
VaR <-quantile(rvec,alpha)
ES <-mean(rvec[rvec<VaR])
round(VaR,6)
round(ES,6)

#sim1 - student-t added 10day
rvec <- rep(0,100000)
for (i in 1:10) {   rvec  <-  rvec+rt.scaled(100000,mean=t.fit$estimate[1],sd=t.fit$estimate[2],df=t.fit$estimate[3]) }
VaR <-quantile(rvec,alpha)
ES <-mean(rvec[rvec<VaR])
round(VaR,6)
round(ES,6)


#sim2 - IID
rvec <- rep(0,100000)
for (i in 1:10) {   rvec <- rvec+ sample(as.vector(logret),100000,replace=TRUE) }
VaR <-quantile(rvec,alpha)
ES <-mean(rvec[rvec<VaR])
round(VaR,6)
round(ES,6)


#sim3- block
rdat <- as.vector(logret) 
rvec <- rep(0,100000) #0 valued vector
posn <- seq(from=1,to=length(rdat)-9,by=1) #position for consecutive from day 1 to last
rpos <- sample(posn,100000,replace=TRUE) #the consecutive day positions
for(i in 1:10) {
  # i-th
  rvec <- rvec+ rdat[rpos]
  rpos <- rpos+1
}
VaR <-quantile(rvec,alpha)
ES <-mean(rvec[rvec<VaR])
round(VaR,6)
round(ES,6)


#auto-correlation
acf(logret)
#auto-correlation absolutely
acf(abs(logret))

#GARCH(1,1)- Normal
uspec <- ugarchspec(  variance.model = list(model = "sGARCH",garchOrder = c(1,1)), mean.model = list(armaOrder = c(0,0), include.mean = TRUE), distribution.model = "normal")
fit.garch <- ugarchfit(spec = uspec, data = logret[,1])

#GARCH(1,1) -t model
uspec <- ugarchspec(variance.model = list(model = "sGARCH",garchOrder = c(1,1)),
                    mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
                    distribution.model = "std")
fit.garch <- ugarchfit(spec = uspec, data = logret[,1])

#the values
round(fit.garch@fit$coef,6)
options(scipen = 999)

#saving
save1 <- cbind( logret[,1], fit.garch@fit$sigma, fit.garch@fit$z ) 
#names(save1) <- c("logret", "s", "z")

#head(save1,3)
#re-autocorr
acf(save1$fit.garch.fit.z)
acf(abs(save1$fit.garch.fit.z))


set.seed(123789)
boot.garch <- ugarchboot(fit.garch,
                         method = c("Partial","Full")[1],#ignore parameter uncertainity
                         sampling = "raw", # draw from standardized residuals
                         n.ahead = 1,   #1 day ahead return 
                         n.bootpred = 100000, #no of simulations
                         solver = "solnp")
rvec <-boot.garch@fseries
VaR <-quantile(rvec,alpha)
ES <-mean(rvec[rvec<VaR])
round(VaR,6)
round(ES,6)


#print the JPM 10K kind of graph
n2016 <- length(logret["1980-01-01/2016-12-31"])

#function needs high processing
roll.garch <- ugarchroll(spec = uspec,
                         data = logret,
                         n.ahead = 1,
                         forecast.length =1,
                         n.start = n2016,
                         refit.every = 1,
                         refit.window = "recursive",
                         calculate.VaR = TRUE,
                         VaR.alpha = alpha,
                         keep.coef = TRUE
                         )

#data(sp500ret)
#head(sp500ret,3)
#head(logret,3)

#additional specification.
spec = ugarchspec(distribution.model = "std")

#fully working, JPM 10K return
mod = ugarchroll(uspec, data = logret, n.ahead = 1, forecast.length =1,
                 n.start = 5000,  refit.every = 100, refit.window = "recursive", 
                 fit.control = list(), solver = "hybrid",
                 calculate.VaR = TRUE, VaR.alpha = 0.01,
                 keep.coef = TRUE)
report(mod, type="VaR", VaR.alpha = 0.01, conf.level = 0.95) 
report(mod, type="fpm")
plot(mod)
#option 4 for VaR