data {
    matrix[96, 2] signatures; // matrix of categories (rows) by signatures (columns)
    int counts[96]; // data = counts per category
}
parameters {
    real<lower=0,upper=1> proportion_somatic;
}
transformed parameters {
    simplex[96] probs;
    vector<lower=0,upper=1>[2] exposures;
    exposures[1] = 1 - proportion_somatic;
    exposures[2] = proportion_somatic;
    probs = signatures * exposures;
    probs = probs / sum(probs);
}
model {
    proportion_somatic ~ uniform(0, 1);
    counts ~ multinomial(probs);
}
