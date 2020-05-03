#include "portfolio_lib.h"

#include <math.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <limits.h>

gsl_rng* initialize_rng_with_seed(unsigned long seed) {

     /* Seed the MC simulation RNG by incrementing a static variable. Ensures
      * That seed will be unique for as many simulations as the maximum size
      * of an unsigned long. Potential bottleneck to parallel performance */
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlxs2);
    gsl_rng_set(rng, seed);
    return rng;
}


void perform_cholesky(gsl_matrix *matrix, const int NUM_ASSETS) {
    printf("HAHAHA!Entering perform_cholesky!\n");

    for(int i = 0; i < 2; i ++) {
        for(int j = 0; j < 2; j++) {
            printf("%f\n", gsl_matrix_get(matrix, i, j));
        }
    }

    gsl_linalg_cholesky_decomp1(matrix);
    printf("HAHAHA!Leaving perform_cholesky!\n");
    /* Make the Cholesky decomposition matrix L lower triangular */
    for (int i = 0; i < NUM_ASSETS; i++) {
        for (int j = i+1; j < NUM_ASSETS; j++) {
            gsl_matrix_set(matrix,i,j,0.0);
        }
    }
}

gsl_vector* corr_norm_rvars(const int NUM_ASSETS, gsl_rng *rng,
        gsl_matrix *cholesky) {

    /* Create a vector of normal random variables */
    gsl_vector *rans = gsl_vector_alloc(NUM_ASSETS);
    for (int i = 0; i < NUM_ASSETS; i++) {
        gsl_vector_set(rans,i,gsl_ran_ugaussian(rng));
    }

    /* Multiply by Cholesky decomposition to generate correlated random
     * variables */
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, cholesky, rans);

    /* Return pointer to the vector of correlated normal random variables */
    return rans;
}

double one_month_portfolio_return(risky_asset assets[],
        const int NUM_ASSETS, gsl_vector *rans) {
    double tot_ret = 0.0;
    /* Add the weighted return for each asset to the total
     * return for the month */
    for (int i = 0; i < NUM_ASSETS; i++) {
        risky_asset curr = assets[i];
        tot_ret += curr.port_weight * one_month_asset_return(curr.mean,
                curr.sigma, gsl_vector_get(rans,i));
    }
    return tot_ret;
}

double one_month_asset_return(double mean, double sigma, double rand_var) {
    /* 1/12 is the change in time for monthly returns */
//    double delta_t = 1.0/12.0;
    double delta_t = 1.0/252.0;
    /* Implementation of the lognormal price process formula */
    return mean * delta_t + sigma * sqrt(delta_t) * rand_var;
}
