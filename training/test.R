xx = seq(0.1,20,0.1)
plot(0,0,xlim=c(0,20), ylim=c(0,1))
for(i in 1:10)
   lines(xx,pZAGA(xx,mu=5, sigma=i/3, nu=0), ylim=c(0,1))  
