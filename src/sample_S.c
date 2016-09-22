/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2005   The R Development Core Team.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* 
 * Relies on the functions for sampling from a truncated normal distribution
 */

#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>

#include "bmea.h"

/* This will sample a single value from a truncated lognormal distribution
 * as the background, then return a new value for the signal 'S' */
double sample_S(double PM_hijk, double log_mean, double log_sd) {
    
    double lim = log(PM_hijk); // Set the limit on the log_scale
    double B = exp(r_righttruncnorm(lim, log_mean, log_sd)); // Call the function r_righttruncnorm
    
    return PM_hijk - B;

}


/* The method using the inverse cdf & which has been abandoned */
/*
double sample_S(double PM_hijk, double p_hijk, double lambda_ik, double delta_ik) {
  
  double new_p = runif(0, p_hijk); // Sample a p_value
  double log_B = qnorm(new_p, lambda_ik, delta_ik, 0, 1); // The inverse CDF
  return lambda_ik;
  //  return PM_hijk - exp(log_B); // The new candidate value for S

}
*/
