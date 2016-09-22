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

#ifndef scanLimit
#define scanLimit	11.090354888959124950675713943331	/* The Limit of Affymetrix scanner resolution i.e. log 2^16 */
#endif

/* Just prototype the functions from the truncnorm.c file */
double ers_a_inf(double a);
double nrs_a_inf(double a);
double ers_a_b(double a, double b);
double nrs_a_b(double a, double b);
double hnrs_a_b(double a, double b);
double urs_a_b(double a, double b);
double r_lefttruncnorm(double a, double mean, double sd);
double r_righttruncnorm(double b, double mean, double sd);
double r_truncnorm(double a, double b, double mean, double sd);

/* Prototype the S sampler as well */
double sample_S(double PM_hijk, double log_mean, double log_sd);

/* Prototype the setDimensions function */

