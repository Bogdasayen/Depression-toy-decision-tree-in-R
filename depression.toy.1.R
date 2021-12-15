############################################################################
# Depression toy example
# Decision tree with three terminal nodes: response, response but relapse, no response
# Assumes a 30 year time horizon (so costs and QALYs of the terminal nodes are over 30 years)
# Performs PSA and generates incremental net benefits, cost effectiveness acceptability
# curves and estimates value of information using nested simulation
############################################################################


############################################################################
# Howard Thom 14-November-2018
# University of Bristol.
# howard.thom@bristol.ac_uk
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
save_samples <- TRUE

# Number of samples
n_samples <- 1000

# Number of treatments (can change but need to specify costs and treatment 
# effects for all treatments)
n_treat <- 3
t_names <- c("No treatment", "CBT", "Antidepressant")

## Utility functions ###############################################
# Logistic link function to convert probabilities to log odds scale
logit <- function(x) {
	return(log(x / (1 - x)))
}
# Inverse of logit to convert log odds to probability scale
expit<-function(x) {
	return(1 / (1 + exp(-x)))
}

# Function to format the results
format_results<-function(x, n_digits = 2) {
	paste(round(mean(x), digits = n_digits),
	      " (", round(quantile(x, probs = 0.025), digits = n_digits), ", ", 
	      round(quantile(x, probs = 0.975), digits = n_digits), ")", sep= "")
}


# Costs for recovery, relapse, and no recovery over 30 year horizon
c_rec <- rnorm(n = n_samples, mean = 1000, sd = 50)
c_rel <- rnorm(n = n_samples, mean=2000, sd = 100)
c_norec<-rnorm(n = n_samples, mean=2500, sd = 125)

# QALYs for recovery, relapse, and no recovery over 30 year horizon
q_rec <- rnorm(n = n_samples, mean = 26, sd = 2)
q_rel <- rnorm(n = n_samples, mean = 23, sd = 3)
q_norec <- rnorm(n = n_samples, mean = 20, sd = 4)

# Set up data structure to store probabilities of recovery 
# and relapse following recovery
p_rel <- p_rec <- matrix(nrow = n_samples, ncol = n_treat)


# Probabilities for no treatment follow beta distributions.
p_rec[, 1] <- rbeta(n = n_samples,  shape1 = 6,  shape2 = 200)
p_rel[, 1] <- rbeta(n = n_samples,  shape1 = 2,  shape2 = 100)

# Log odds ratios for comparator treatments relative to placebo
# 1 = CBT,  2 = Antidepressant
# CBT has lower recovery but lower relapse probability than antidepressant
# Note that this is not based on any data,  and is only for illustrative purposes

# Log odds ratios of recovery and relapse
# These were estimated using a BUGS network meta-analysis applied to artificial data
lor_rec <- mvrnorm(n = n_samples, mu = c(0.99, 1.33), 
                   Sigma = matrix(c(0.22, 0.15, 0.15, 0.20), nrow = 2))
lor_rel <- mvrnorm(n = n_samples, mu = c(-1.48, -0.40), 
                   Sigma = matrix(c(0.14, 0.05, 0.05, 0.11), nrow = 2))

# Cost of treatment,  so no treatment is free,  CBT is expensive,  antidepressants are cheap
# CBT is approximate cost for 10 sessions at £30 per session. 
c_treat <- t(matrix(rep(c(0, 300, 30), n_samples), ncol = n_samples, nrow = 3))

# Willingness to pay thresholds 
# Set it to a vector so we can draw the cost-effectiveness acceptability curve (CEAC)
lambdas <- c(1:50) * 1000
# Lamdba target is the key threshold for decision making
lambda.target <- 20000

# Build matrices to stroke absolute / incremental costs,  effects,  and net benefits
incremental_costs <- incremental_effects <- incremental_nb <- costs <- effects <- net_benefit <- matrix(nrow = n_samples, ncol = n_treat)
# Name the columns after the treatments
colnames(incremental_costs) <- colnames(incremental_effects) <- colnames(incremental_nb) <- colnames(p_rec) <- colnames(p_rel) <- colnames(effects) <- colnames(costs) <- colnames(net_benefit) <- t_names

# Use the absolute probabilities of recovery and relapse for no treatment
# and the log odds ratios for CBT and antidepressants
for(i in 2:3)
{
  p_rec[, i] <- expit(logit(p_rec[, 1]) + lor_rec[, i - 1])
  p_rel[, i] <- expit(logit(p_rel[, 1]) + lor_rel[, i - 1])
}

# This can be vectorised as below,  but only a loop over two treatments so
# speed advantage is limited,  while clarity is lost. Note that the code
# is already vectorised over the PSA samples.
# p_rec[, c(2:n_treat)] <- expit(logit(p_rec[, 1]) + lor_rec[, c(2:n_treat) - 1])
# p_rel[, c(2:n_treat)] <- expit(logit(p_rel[, 1]) + lor_rel[, c(2:n_treat) - 1])


# The following two lines are the entire decision tree calculation
# Add extra probabilities for additional branches to the model
# This can be deterministic or probabilistic,  the code is the same.
effects <- p_rec * (1-p_rel) * q_rec + p_rec * p_rel * q_rel + (1-p_rec) * q_norec
costs <- c_treat + p_rec * (1-p_rel) * c_rec + p_rec * p_rel * c_rel + (1-p_rec) * c_norec

# Now calculate the net benefit at "lambda target"
net_benefit <- lambda.target * effects - costs

# Incremental costs,  effects,  and net benefits
incremental_costs <- costs-costs[, 1]
incremental_effects <- effects-effects[, 1]
incremental_nb <- net_benefit - net_benefit[, 1]

# Can use the colMeans() function to get a quick look at the
# average results_ These are point estimates.
colMeans(costs)
colMeans(effects)
colMeans(net_benefit)

# If you want to export results to SAVI / BCEA to estimate EVPPI or do
# other anlaysis,  use the following code.
# Build a general input_parameters matrix
input_parameters <- t(rbind(t(p_rec), t(p_rel), t(lor_rec), t(lor_rel), q_rec, q_rel, q_norec, c_rec, c_rel, c_norec))
colnames(input_parameters)[1:6] <- c(paste("p_rec", colnames(p_rec)), paste("p_rel", colnames(p_rel)))
colnames(input_parameters)[7:10] <- c(paste("lor_rec", colnames(p_rec)[2:3]), paste("lor_rel", colnames(p_rel)[2:3]))

# Export data for SAVI / BCEA
#save(costs, effects, input_parameters, file = paste(baseline.directory, "savi_bcea.data.", n_samples, ".rda", sep = ""))
#write.csv(costs, file = paste(baseline.directory, "costs.", n_samples, ".csv", sep = ""))
#write.csv(effects, file = paste(baseline.directory, "effects.", n_samples, ".csv", sep = ""))
#write.csv(input_parameters, file = paste(baseline.directory, "input_parameters.", n_samples, ".csv", sep = ""))


# Build a matrix to store the results
results_matrix <- matrix(NA,  nrow = 6, ncol = n_treat)
rownames(results_matrix) <- c("Total costs", "Total QALYs", "Incremental costs", "Incremental QALYs", "Net Benefit", "Incremental NB")
colnames(results_matrix) <- t_names
for(i_treat in 1:n_treat)
{
  results_matrix["Total costs", i_treat] <- format_results(x = costs[, i_treat])
  results_matrix["Total QALYs", i_treat] <- format_results(x = effects[, i_treat])
  results_matrix["Incremental costs", i_treat] <- format_results(x = incremental_costs[, i_treat])
  results_matrix["Incremental QALYs", i_treat] <- format_results(x = incremental_effects[, i_treat])
  results_matrix["Net Benefit", i_treat] <- format_results(x = net_benefit[, i_treat])
  results_matrix["Incremental NB", i_treat] <- format_results(x = incremental_nb[, i_treat])
}

# Export as a csv
write.csv(results_matrix, file = "depression_results.csv")
# Or as an Excel file
write.xlsx(results_matrix, file = "depression_results.xlsx", sheetName = "CEA results")

# Vector of optimal treatment at £20, 000 for each PSA
which_max_nb <- apply(net_benefit, c(1), which.max)


# Calculate  CEAC
net_benefits <- array(NA,  dim = c(length(lambdas), dim(effects)[1], dim(effects)[2]))
ceac <- matrix(nrow = length(lambdas), ncol = n_treat)
for(i_lambda in 1:length(lambdas))
{
  # net benefit for each threshold i_lambda
  net_benefits[i_lambda, , ] <- lambdas[i_lambda] * effects-costs
  which_max_nb <- apply(net_benefits[i_lambda, , ], c(1), which.max)
  for(i_treat in 1:n_treat)
  {
    # Probability i_treat is optimal at i_lambda
    ceac[i_lambda, i_treat] <- mean(which_max_nb == i_treat)
  }
}
# Use the following line if you want to export to a jpeg file (or use pdf() for same result)
# Can't use both simultaneously and be sure to use dev.off() to close connection to file
#jpeg(file = paste(baseline.directory, "ceac_depression.toy.", n_samples, ".jpg", sep = ""))
#pdf(file = paste(baseline.directory, "ceac_depression.toy.", n_samples, ".pdf", sep = ""))

# Plot the CEAC
plot(c(0, 0), col = 0, xlim = c(0, max(lambdas)), ylim = c(0, 1.05), main = "Cost-effectiveness acceptability curve", xlab = "Willingness-to-pay (£)", ylab = "Probability most cost-effective")
for(i_treat in 1:n_treat)
{
  lines(lambdas,  ceac[, i_treat], lty = i_treat, lwd = 3, col = i_treat)
}
legend("topleft", legend = t_names, lty = c(1:n_treat), lwd = 3, col = c(1:n_treat))

# Call this to close access to the jpeg() or pdf()
#dev.off()

## Value of information analysis #############################################
# Net benefit function so this is not expensive to estimate
# Can compare with BCEA / SAVI estimates
# Calculate the total EVPI
evpi <- mean(apply(net_benefit, c(1), max))-mean(net_benefit[, which.max(colMeans(net_benefit))])


# Need nested simulation for EVPPI
# Epidimiological focal parameters p_rec,  p_rel
effects_inner <- p_rec * (1-p_rel) * mean(q_rec) + p_rec * p_rel * mean(q_rel) + (1-p_rec) * mean(q_norec)
costs_inner <- c_treat + p_rec * (1-p_rel) * mean(c_rec) + p_rec * p_rel * mean(c_rel) + (1-p_rec) * mean(c_norec)
net_benefit_inner <- lambda.target * effects_inner-costs_inner
evppi_epi <- mean(apply(net_benefit_inner, c(1), max))-mean(net_benefit[, which.max(colMeans(net_benefit))])


# Need nested simulation for EVPPI
# Cost-effectiveness focal parameters q_rec,  q_rel,  q_norec,  c_rec,  c_rel,  c_norec
effects_inner <- t(outer(colMeans(p_rec) * (1-colMeans(p_rel)), q_rec) + outer(colMeans(p_rec) * colMeans(p_rel), q_rel) + outer((1-colMeans(p_rec)), q_norec))
costs_inner <- c_treat + t(outer(colMeans(p_rec) * (1-colMeans(p_rel)), c_rec) + outer(colMeans(p_rec) * colMeans(p_rel), c_rel) + outer((1-colMeans(p_rec)), c_norec))
net_benefit_inner <- lambda.target * effects_inner-costs_inner
evppi_ce <- mean(apply(net_benefit_inner, c(1), max))-mean(net_benefit[, which.max(colMeans(net_benefit))])

# Need nested simulation for EVPPI
# CBT focal parameters p_rec[, 2],  p_rel[, 2]
n_samples_inner <- 100 

rec_mean <- colMeans(lor_rec); rec_sd <- sqrt(c(var(lor_rec)[1, 1], var(lor_rec)[2, 2]))
rec_corr <- cor(lor_rec[, 1], lor_rec[, 2])
rel.mean <- colMeans(lor_rel); rel.sd <- (c(var(lor_rel)[1, 1], var(lor_rel)[2, 2]))
rel.corr <- cor(lor_rel[, 1], lor_rel[, 2])

# Antidepressant lors depend on CBT lor so need conditional sampling
lor_rec_inner <- lor_rel_inner <- array(NA,  dim = c(n_samples_inner, 2))
p_rec_inner <- p_rel_inner <- array(NA,  dim = c(n_samples_inner,  3))
effects_inner <- costs_inner <- matrix(NA, nrow = n_samples, ncol = 3)

# Loop over the outer samples

for(i in 1:n_samples)
{	
  # Take inner samples for each outer sample (not very memory efficient)
  # Probabilities for no treatment
  p_rec_inner[, 1] <- rbeta(n = n_samples_inner,  shape1 = 6,  shape2 = 200)
  p_rel_inner[, 1] <- rbeta(n = n_samples_inner,  shape1 = 2,  shape2 = 100)
  
  lor_rec_inner[, 1] <- lor_rec[i, 1]
  lor_rel_inner[, 1] <- lor_rel[i, 1]
  lor_rec_inner[, 2] <- rnorm(n = n_samples_inner, mean = rec_mean[2] + (rec_sd[2] / rec_sd[1]) * rec_corr * (lor_rec[i, 1]-rec_mean[1]), sd = sqrt((1-rec_corr^2) * rec_sd[2]^2))
  lor_rel_inner[, 2] <- rnorm(n = n_samples_inner, mean = rel.mean[2] + (rel.sd[2] / rel.sd[1]) * rel.corr * (lor_rel[i, 1]-rel.mean[1]), sd = sqrt((1-rel.corr^2) * rel.sd[2]^2))
  for(j in 2:3)
  {
    p_rec_inner[, j] <- expit(logit(p_rec_inner[, 1]) + lor_rec_inner[, j-1])
    p_rel_inner[, j] <- expit(logit(p_rel_inner[, 1]) + lor_rel_inner[, j-1])
  }
  
  # Calculate the inner expectation for each outer sample
  for(j in 1:3)
  {
    effects_inner[i, j] <- mean(p_rec_inner[, j]) * (1-mean(p_rel_inner[, j])) * mean(q_rec) + mean(p_rec_inner[, j]) * mean(p_rel_inner[, j]) * mean(q_rel) + (1-mean(p_rec_inner[, j])) * mean(q_norec)
    costs_inner[i, j] <- c_treat[j] + mean(p_rec_inner[, j]) * (1-mean(p_rel_inner[, j])) * mean(c_rec) + mean(p_rec_inner[, j]) * mean(p_rel_inner[, j]) * mean(c_rel) + (1-mean(p_rec_inner[, j])) * mean(c_norec)
  }
}

# And finally calculate the EVPPI of the CBT
net_benefit_inner <- lambda.target * effects_inner-costs_inner
evppi_cbt <- mean(apply(net_benefit_inner, c(1), max))-mean(net_benefit[, which.max(colMeans(net_benefit))])




