functions {
  // function to calculate the pdf of 
  // a special prior distribution for the beta_g
  // sb is the (positive) scale parameter of a normal distribution
  // c is the (positive) cutoff for minimum log fc - any absolute value;
  //  smaller than the cutoff will have a density of 0
  real cbp_lpdf(real x, real sb, real c) {
    real lnorm = log(1 - (normal_cdf(c | 0, sb) - normal_cdf(-c | 0, sb)));
    real out;
    if (abs(x) > c) {
      out = normal_lpdf(x | 0, sb) - lnorm;
    } else if (abs(x) <= c) {
      out = 0;
    }
    return(out);
  }
}