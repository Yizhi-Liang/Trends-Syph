data {
  
  int<lower=1> n; // Total number of observations
  int<lower=1> n_year; // Number of unique years
  int<lower=1, upper=n_year> year[n]; // Array of year for each observation
  int<lower=1> n_race; // Number of unique races
  int<lower=1, upper=n_race> race[n]; // Array of race for each observation
  int<lower=0> Positivities[n]; // Number of positivities for each observation
  int<lower=0> Pop[n]; // Population for each observation
  int<lower=0> Pop_adj[n]; // Adjusted Population size
  real<lower=0> Test_alpha[n]; // Prior alpha for proportion tested for each observation
  real<lower=0> Test_beta[n]; // Prior beta for proportion tested for each observation
  real<lower=0> Sens_alpha; // Prior alpha for sensitivity for each observation
  real<lower=0> Sens_beta; // Prior beta for sensitivity for each observation
  real<lower=0> Spec_alpha; // Prior alpha for specificity for each observation
  real<lower=0> Spec_beta; // Prior beta for specificity for each observation
  real<lower=0> P_alpha; // Prior alpha for prevalence for each observation
  real<lower=0> P_beta; // Prior beta for prevalence for each observation
  
}

parameters {
  
  real<lower=0, upper=1> P[n_year, n_race]; // True prevalence for each year and race
  real<lower=0, upper=1> Sens;              // Sensitivity
  real<lower=0, upper=1> Spec;              // Specificity
  real<lower=0, upper=1> Test[n_year, n_race];         // Proportion tested for each race
  
}

transformed parameters {
  
  real<lower=0> theta[n_year, n_race];  // Expected number of positive tests for each year and race
  
  for (i in 1:n) {
    theta[year[i], race[i]] = ((P[year[i], race[i]] * Sens) + (1 - P[year[i], race[i]]) * (1 - Spec)) * 
    Test[year[i], race[i]];
  }
  
}

model {
  
  // Priors
  
  for (i in 1:n) {
    P[year[i], race[i]] ~ beta(P_alpha, P_beta); // Prior for prevalence (P)

    Test[year[i], race[i]] ~ beta(Test_alpha[i], Test_beta[i]); // Prior for proportion tested (Test）
  }
  
  Sens ~ beta(Sens_alpha, Sens_beta); // Prior for sensitivity (Sens)
  Spec ~ beta(Spec_alpha, Spec_beta); // Prior for specificity (Spec)
  
  // Likelihood of the observed data given the model
  for (i in 1:n) {
    Positivities[i] ~ binomial(Pop_adj[i], theta[year[i], race[i]]);
  }
  
}

generated quantities {

  int<lower=0> pred_Positivities[n_year, n_race]; // Predicted number of positive tests for each year and race

  //  predictions
  for (i in 1:n) {
    pred_Positivities[year[i], race[i]] = binomial_rng(Pop_adj[i], theta[year[i], race[i]]);
  }

}
