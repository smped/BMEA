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
 *Contains the functions for sampling from a truncated normal distribution
 */

/* 
 * Based almost entirely on those written for the truncnorm package by
 * Björn Bornkamp   <bornkamp@statistik.tu-dortmund.de>
 * Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
 */

#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>

/* The required variables for the truncated normal sampling are: */
const double t1 = 0.15;
const double t2 = 2.18;
const double t3 = 0.725;
const double t4 = 0.45;

/* Exponential rejection sampling (a,inf) */
double ers_a_inf(double a) {
    const double ainv = 1.0 / a;
    double x,rho;
    do{
	x = rexp(ainv) + a; /* rexp works with 1/lambda */
	rho = exp(-0.5 * pow((x - a), 2));
    } while (runif(0, 1) > rho);
    return x;
}

/* Normal rejection sampling (a,inf) */
double nrs_a_inf(double a){
  double x = a - 1.0;
  while(x < a){
    x = rnorm(0, 1);
  }
  return x;
}

/* Exponential rejection sampling (a,b) */
double ers_a_b(double a, double b) {
    const double ainv = 1.0 / a;
    double x, rho;
    do{
	x = rexp(ainv) + a; /* rexp works with 1/lambda */
	rho = exp(-0.5 * pow((x-a), 2));
    } while (runif(0, 1) > rho || x > b);
    return x;
}

/* Normal rejection sampling (a,b) */
double nrs_a_b(double a, double b){
  double x = a - 1.0;
  while(x < a || x > b){
    x = rnorm(0, 1);
  }
  return x;
}

/* Half-normal rejection sampling */
double hnrs_a_b(double a, double b){
  double x = a - 1.0;
  while(x < a || x > b) {
    x = rnorm(0, 1);
    x = fabs(x);
  }
  return x;
}

/* Uniform rejection sampling */
double urs_a_b(double a, double b){
    const double phi_a = dnorm(a, 0.0, 1.0, FALSE);
    double x = 0.0, u = 0.0;
    if (a < 0 && b > 0) {
	do{
	    x = runif(a, b);
	} while (runif(0, 1) * M_1_SQRT_2PI > dnorm(x,0,1,0));
    } else {
	do{
	    x = runif(a, b);
	    u = runif(0, 1);
	} while (runif(0, 1) * phi_a > dnorm(x, 0.0, 1.0, FALSE));
    }
    return x;
}


/* Previously this was refered to as type 1 sampling: */
double r_lefttruncnorm(double a, double mean, double sd) {

  const double alpha = (a - mean) / sd;
  if (alpha < t4) { // If alpha<t4, use normal rejection sampling
    return mean + sd * nrs_a_inf(alpha);
  } else { // Otherwise use exponential rejection sampling
    return mean + sd * ers_a_inf(alpha);
  }
  
  // Use the inverse cdf method, which is slower...
  /*
  double lim = pnorm(a, mean, sd, 1, 0); // Find the p-value of the lower limit
  double p = runif(lim, 1); // Sample from U(lim, 1)
  double x = qnorm(p, mean, sd, 1, 0); // Use the inverse cdf transformation for the random value
  return x;
  */

}

double r_righttruncnorm(double b, double mean, double sd) {
  const double beta = (b - mean) / sd;
  //Exploit symmetry: 
  return mean - sd * r_lefttruncnorm(-beta, 0.0, 1.0); // Call the function r_lefttruncnorm

  // Don't use the inverse cdf method, which is slower
  /*
  double lim = pnorm(b, mean, sd, 1, 0); // Find the p-value of the upper limit
  double p = runif(0, lim); // Sample from U(0, lim)
  double x = qnorm(p, mean, sd, 1, 0); // Use the inverse cdf transformation for the random value
  return x;
  */
}

double r_truncnorm(double a, double b, double mean, double sd) {
  
  const double alpha = (a - mean) / sd;
  const double beta = (b - mean) / sd;
  const double phi_a = dnorm(a, 0.0, 1.0, FALSE);
  const double phi_b = dnorm(b, 0.0, 1.0, FALSE);
  if (beta <= alpha) {
    return NA_REAL;
  } else if (alpha < 0 && 0 < beta) {
    if (phi_a < t1 && phi_b < t1) {
      return mean + sd * nrs_a_b(alpha, beta);
    } else {
      return mean + sd * urs_a_b(alpha, beta);
    }
  } else if (alpha >= 0) {
    if (phi_a / phi_b <= t2) {
      return mean + sd * urs_a_b(alpha, beta);
    } else {
      if (alpha < t3) {
	return mean + sd * hnrs_a_b(alpha, beta);
      } else {
	return mean + sd * ers_a_b(alpha, beta);
      }
    }
  } else { /* beta <= 0 */
    if (phi_a / phi_b <= t2) {
      return mean - sd * urs_a_b(-beta, -alpha);
    } else {
      if (alpha < t3) {
	return mean - sd * hnrs_a_b(-beta, -alpha);
      } else {
	return mean - sd * ers_a_b(-beta, -alpha);
      }
    }
  }

  /* The slower inverse cdf method
  if (b <= a) { // Check the limits are valis
    return NA_REAL;
  }
  else{
    double lowLim = pnorm(a, mean, sd, 1, 0); // Find the p-value of the lower limit
    double upLim = pnorm(b, mean, sd, 1, 0);
    double p = runif(lowLim, upLim); // Sample from U(lowLim, upLim)
    double x = qnorm(p, mean, sd, 1, 0); // Use the inverse cdf transformation for the random value
    return x;
  }
  */
  
}
