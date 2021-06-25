############################################################################
# Depression toy example
# Decision tree with three terminal nodes: response, response but relapse, no response
# Assumes a 30 year time horizon (so costs and QALYs of the terminal nodes are over 30 years)
# Performs PSA and generates incremental net benefits, cost effectiveness acceptability
# curves and estimates value of information using nested simulation
############################################################################

# Modification to test GitHub desktop

############################################################################
# Howard Thom 14-November-2018
# University of Bristol.
# howard.thom@bristol.ac.uk
############################################################################

############################################################################
# The HTMR Collaboration and innovation in Difficult and Complex randomised 
# controlled Trials In Invasive procedures (ConDuCT-II) hub provided support for this research.
# This study was supported by the NIHR Biomedical Research Centre at the 
# University Hospitals Bristol NHS Foundation Trust and the University of 
# Bristol. The views expressed in this presentation are those of the author(s)
# and not necessarily those of the NHS, the National Institute for Health Research or the Department of Health.
############################################################################


# Set to desired baseline directory
baseline.directory<-"C:/Users/Howard/Documents/Bristol/R for CEA/ISPOR workshop Decision models with R 2018/code examples"
setwd(baseline.directory)


# For mvrnorm function
library(MASS)
# For reading/writing excel files
library(xlsx)

# Save inputs and outputs to file?
save.samples<-TRUE

# Number of samples
n.samples<-10000

# Number of treatments (can change but need to specify costs and treatment 
# effects for all treatments)
n.treat<-3
t.names<-c("No treatment","CBT","Antidepressant")

## Utility functions ###############################################
# Logistic link function to convert probabilities to log odds scale
logit<-function(x)
{
	return(log(x/(1-x)))
}
# Inverse of logit to convert log odds to probability scale
expit<-function(x)
{
	return(1/(1+exp(-x)))
}
# Function to format the results
format.results<-function(x,n.digits=2)
{
	paste(round(mean(x),digits=n.digits)," (",round(quantile(x,probs=0.025),digits=n.digits),", ",round(quantile(x,probs=0.975),digits=n.digits),")",sep="")
}


# Costs for recovery, relapse, and no recovery over 30 year horizon
c.rec<-rnorm(n=n.samples, mean=1000, sd=50)
c.rel<-rnorm(n=n.samples, mean=2000, sd=100)
c.norec<-rnorm(n=n.samples, mean=2500, sd=125)

# QALYs for recovery, relapse, and no recovery over 30 year horizon
q.rec<-rnorm(n=n.samples, mean=26, sd=2)
q.rel<-rnorm(n=n.samples, mean=23, sd=3)
q.norec<-rnorm(n=n.samples, mean=20, sd=4)

# Set up data structure to store probabilities of recovery 
# and relapse following recovery
p.rel<-p.rec<-matrix(nrow=n.samples, ncol=n.treat)

# Probabilities for no treatment follow beta distributions.
p.rec[,1]<-rbeta(n=n.samples, shape1=6, shape2=200)
p.rel[,1]<-rbeta(n=n.samples, shape1=2, shape2=100)

# Log odds ratios for comparator treatments relative to placebo
# 1=CBT, 2=Antidepressant
# CBT has lower recovery but lower relapse probability than antidepressant
# Note that this is not based on any data, and is only for illustrative purposes

# Log odds ratios of recovery and relapse
# These were estimated using a BUGS network meta-analysis applied to artificial data
lor.rec<-mvrnorm(n=n.samples,mu=c(0.99,1.33),Sigma=matrix(c(0.22,0.15,0.15,0.20),nrow=2))
lor.rel<-mvrnorm(n=n.samples,mu=c(-1.48,-0.40),Sigma=matrix(c(0.14,0.05,0.05,0.11),nrow=2))

# Cost of treatment, so no treatment is free, CBT is expensive, antidepressants are cheap
# CBT is approximate cost for 10 sessions at £30 per session. 
c.treat<-t(matrix(rep(c(0,300,30),n.samples),ncol=n.samples,nrow=3))

# Willingness to pay thresholds 
# Set it to a vector so we can draw the cost-effectiveness acceptability curve (CEAC)
lambdas<-c(1:50)*1000
# Lamdba target is the key threshold for decision making
lambda.target<-20000

# Build matrices to stroke absolute/incremental costs, effects, and net benefits
incremental.costs<-incremental.effects<-incremental.nb<-costs<-effects<-net.benefit<-matrix(nrow=n.samples,ncol=n.treat)
# Name the columns after the treatments
colnames(incremental.costs)<-colnames(incremental.effects)<-colnames(incremental.nb)<-colnames(p.rec)<-colnames(p.rel)<-colnames(effects)<-colnames(costs)<-colnames(net.benefit)<-t.names

# Use the absolute probabilities of recovery and relapse for no treatment
# and the log odds ratios for CBT and antidepressants
for(i in 2:3)
{
	p.rec[,i]<-expit(logit(p.rec[,1])+lor.rec[,i-1])
	p.rel[,i]<-expit(logit(p.rel[,1])+lor.rel[,i-1])
}

# This can be vectorised as below, but only a loop over two treatments so
# speed advantage is limited, while clarity is lost. Note that the code
# is already vectorised over the PSA samples.
# p.rec[,c(2:n.treat)]<-expit(logit(p.rec[,1])+lor.rec[,c(2:n.treat)-1])
# p.rel[,c(2:n.treat)]<-expit(logit(p.rel[,1])+lor.rel[,c(2:n.treat)-1])


# The following two lines are the entire decision tree calculation
# Add extra probabilities for additional branches to the model
# This can be deterministic or probabilistic, the code is the same.
effects<-p.rec*(1-p.rel)*q.rec+p.rec*p.rel*q.rel+(1-p.rec)*q.norec
costs<-c.treat+p.rec*(1-p.rel)*c.rec+p.rec*p.rel*c.rel+(1-p.rec)*c.norec

# Now calculate the net benefit at "lambda target"
net.benefit<-lambda.target*effects-costs

# Incremental costs, effects, and net benefits
incremental.costs<-costs-costs[,1]
incremental.effects<-effects-effects[,1]
incremental.nb<-net.benefit-net.benefit[,1]

# Can use the colMeans() function to get a quick look at the
# average results. These are point estimates.
colMeans(costs)
colMeans(effects)
colMeans(net.benefit)

# If you want to export results to SAVI/BCEA to estimate EVPPI or do
# other anlaysis, use the following code.
# Build a general input.parameters matrix
input.parameters<-t(rbind(t(p.rec),t(p.rel),t(lor.rec),t(lor.rel),q.rec,q.rel,q.norec,c.rec,c.rel,c.norec))
colnames(input.parameters)[1:6]<-c(paste("p.rec",colnames(p.rec)),paste("p.rel",colnames(p.rel)))
colnames(input.parameters)[7:10]<-c(paste("lor.rec",colnames(p.rec)[2:3]),paste("lor.rel",colnames(p.rel)[2:3]))

# Export data for SAVI/BCEA
#save(costs,effects,input.parameters,file=paste(baseline.directory,"savi.bcea.data.",n.samples,".rda",sep=""))
#write.csv(costs,file=paste(baseline.directory,"costs.",n.samples,".csv",sep=""))
#write.csv(effects,file=paste(baseline.directory,"effects.",n.samples,".csv",sep=""))
#write.csv(input.parameters,file=paste(baseline.directory,"input.parameters.",n.samples,".csv",sep=""))


# Build a matrix to store the results
results.matrix<-matrix(NA, nrow=6,ncol=n.treat)
rownames(results.matrix)<-c("Total costs","Total QALYs","Incremental costs","Incremental QALYs","Net Benefit","Incremental NB")
colnames(results.matrix)<-t.names
for(i.treat in 1:n.treat)
{
	results.matrix["Total costs",i.treat]<-format.results(x=costs[,i.treat])
	results.matrix["Total QALYs",i.treat]<-format.results(x=effects[,i.treat])
	results.matrix["Incremental costs",i.treat]<-format.results(x=incremental.costs[,i.treat])
	results.matrix["Incremental QALYs",i.treat]<-format.results(x=incremental.effects[,i.treat])
	results.matrix["Net Benefit",i.treat]<-format.results(x=net.benefit[,i.treat])
	results.matrix["Incremental NB",i.treat]<-format.results(x=incremental.nb[,i.treat])
}

# Export as a csv
write.csv(results.matrix,file="depression.results.csv")
# Or as an Excel file
write.xlsx(results.matrix,file="depression.results.xlsx",sheetName="CEA results")

# Vector of optimal treatment at £20,000 for each PSA
which.max.nb<-apply(net.benefit,c(1),which.max)


# Calculate  CEAC
net.benefits<-array(NA, dim=c(length(lambdas),dim(effects)[1],dim(effects)[2]))
ceac<-matrix(nrow=length(lambdas),ncol=n.treat)
for(i.lambda in 1:length(lambdas))
{
	# net benefit for each threshold i.lambda
	net.benefits[i.lambda,,]<-lambdas[i.lambda]*effects-costs
	which.max.nb<-apply(net.benefits[i.lambda,,],c(1),which.max)
	for(i.treat in 1:n.treat)
	{
		# Probability i.treat is optimal at i.lambda
		ceac[i.lambda,i.treat]<-mean(which.max.nb==i.treat)
	}
}
# Use the following line if you want to export to a jpeg file (or use pdf() for same result)
# Can't use both simultaneously and be sure to use dev.off() to close connection to file
#jpeg(file=paste(baseline.directory,"ceac.depression.toy.",n.samples,".jpg",sep=""))
#pdf(file=paste(baseline.directory,"ceac.depression.toy.",n.samples,".pdf",sep=""))

# Plot the CEAC
plot(c(0,0),col=0,xlim=c(0,max(lambdas)),ylim=c(0,1.05),main="Cost-effectiveness acceptability curve",xlab="Willingness-to-pay (£)",ylab="Probability most cost-effective")
for(i.treat in 1:n.treat)
{
	lines(lambdas, ceac[,i.treat],lty=i.treat,lwd=3,col=i.treat)
}
legend("topleft",legend=t.names,lty=c(1:n.treat),lwd=3,col=c(1:n.treat))

# Call this to close access to the jpeg() or pdf()
#dev.off()

## Value of information analysis #############################################
# Net benefit function so this is not expensive to estimate
# Can compare with BCEA/SAVI estimates
# Calculate the total EVPI
evpi<-mean(apply(net.benefit,c(1),max))-mean(net.benefit[,which.max(colMeans(net.benefit))])


# Need nested simulation for EVPPI
# Epidimiological focal parameters p.rec, p.rel
effects.inner<-p.rec*(1-p.rel)*mean(q.rec)+p.rec*p.rel*mean(q.rel)+(1-p.rec)*mean(q.norec)
costs.inner<-c.treat+p.rec*(1-p.rel)*mean(c.rec)+p.rec*p.rel*mean(c.rel)+(1-p.rec)*mean(c.norec)
net.benefit.inner<-lambda.target*effects.inner-costs.inner
evppi.epi<-mean(apply(net.benefit.inner,c(1),max))-mean(net.benefit[,which.max(colMeans(net.benefit))])


# Need nested simulation for EVPPI
# Cost-effectiveness focal parameters q.rec, q.rel, q.norec, c.rec, c.rel, c.norec
effects.inner<-t(outer(colMeans(p.rec)*(1-colMeans(p.rel)),q.rec)+outer(colMeans(p.rec)*colMeans(p.rel),q.rel)+outer((1-colMeans(p.rec)),q.norec))
costs.inner<-c.treat+t(outer(colMeans(p.rec)*(1-colMeans(p.rel)),c.rec)+outer(colMeans(p.rec)*colMeans(p.rel),c.rel)+outer((1-colMeans(p.rec)),c.norec))
net.benefit.inner<-lambda.target*effects.inner-costs.inner
evppi.ce<-mean(apply(net.benefit.inner,c(1),max))-mean(net.benefit[,which.max(colMeans(net.benefit))])

# Need nested simulation for EVPPI
# CBT focal parameters p.rec[,2], p.rel[,2]
n.samples.inner<-100 

rec.mean<-colMeans(lor.rec); rec.sd<-sqrt(c(var(lor.rec)[1,1],var(lor.rec)[2,2]))
rec.corr<-cor(lor.rec[,1],lor.rec[,2])
rel.mean<-colMeans(lor.rel); rel.sd<-(c(var(lor.rel)[1,1],var(lor.rel)[2,2]))
rel.corr<-cor(lor.rel[,1],lor.rel[,2])

# Antidepressant lors depend on CBT lor so need conditional sampling
lor.rec.inner<-lor.rel.inner<-array(NA, dim=c(n.samples.inner,2))
p.rec.inner<-p.rel.inner<-array(NA, dim=c(n.samples.inner, 3))
effects.inner<-costs.inner<-matrix(NA,nrow=n.samples,ncol=3)

# Loop over the outer samples

for(i in 1:n.samples)
{	
	# Take inner samples for each outer sample (not very memory efficient)
	# Probabilities for no treatment
	p.rec.inner[,1]<-rbeta(n=n.samples.inner, shape1=6, shape2=200)
	p.rel.inner[,1]<-rbeta(n=n.samples.inner, shape1=2, shape2=100)

	lor.rec.inner[,1]<-lor.rec[i,1]
	lor.rel.inner[,1]<-lor.rel[i,1]
	lor.rec.inner[,2]<-rnorm(n=n.samples.inner,mean=rec.mean[2]+(rec.sd[2]/rec.sd[1])*rec.corr*(lor.rec[i,1]-rec.mean[1]),sd=sqrt((1-rec.corr^2)*rec.sd[2]^2))
	lor.rel.inner[,2]<-rnorm(n=n.samples.inner,mean=rel.mean[2]+(rel.sd[2]/rel.sd[1])*rel.corr*(lor.rel[i,1]-rel.mean[1]),sd=sqrt((1-rel.corr^2)*rel.sd[2]^2))
	for(j in 2:3)
	{
		p.rec.inner[,j]<-expit(logit(p.rec.inner[,1])+lor.rec.inner[,j-1])
		p.rel.inner[,j]<-expit(logit(p.rel.inner[,1])+lor.rel.inner[,j-1])
	}

	# Calculate the inner expectation for each outer sample
	for(j in 1:3)
	{
		effects.inner[i,j]<-mean(p.rec.inner[,j])*(1-mean(p.rel.inner[,j]))*mean(q.rec)+mean(p.rec.inner[,j])*mean(p.rel.inner[,j])*mean(q.rel)+(1-mean(p.rec.inner[,j]))*mean(q.norec)
		costs.inner[i,j]<-c.treat[j]+mean(p.rec.inner[,j])*(1-mean(p.rel.inner[,j]))*mean(c.rec)+mean(p.rec.inner[,j])*mean(p.rel.inner[,j])*mean(c.rel)+(1-mean(p.rec.inner[,j]))*mean(c.norec)
	}
}

# And finally calculate the EVPPI of the CBT
net.benefit.inner<-lambda.target*effects.inner-costs.inner
evppi.cbt<-mean(apply(net.benefit.inner,c(1),max))-mean(net.benefit[,which.max(colMeans(net.benefit))])




