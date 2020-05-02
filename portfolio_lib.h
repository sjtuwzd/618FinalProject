#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

#ifndef PORTFOLIO_LIB_H
#define PORTFOLIO_LIB_H
/* Stores statistics information on the annual distribution of returns 
 * of a single risky asset along with its ticker and
 * weight within the portfolio */
struct {
    char *ticker;
    double mean;
    double sigma; /* The standard deviation of returns */
    double port_weight;
} typedef risky_asset;

/* Create a new random number generator using the specified seed
 *
 * \param seed is the seed to initialize the new RNG with
 * \return a new random number generator initialized with the given seed
 *
 * Notes:
 * The memory of the returned RNG pointer is allocated inside this function.
 * The caller is responsible for freeing that memory
 */
gsl_rng* initialize_rng_with_seed(unsigned long seed);

/* Evaluate the Cholesky decomposition on a variance-covariance matrix
 *
 * \param matrix is the variance-covariance matrix to perform the Cholesky
 *  decomposition on; the contents of this matrix will be replaced with that
 *  of the Cholesky decomposition matrix
 * \param NUM_ASSETS is the number of assets represented in the
 *  variance-covariance matrix
 */
void perform_cholesky(gsl_matrix *matrix, const int NUM_ASSETS);

/* Returns a pointer to a vector containing correlated, normally distributed 
 * random variables given the Cholesky decomposition
 *
 * \param NUM_ASSETS the number of assets to be generated random variables for
 * \param rng a pointer to the random number generator that is sampled to
 *  generate the random variables
 * \param cholesky the Cholesky decomposition matrix of the assets
 * \return a pointer to a vector of NUM_ASSETS correlated random variables
 *
 * Notes:
 * The returned vector pointer is allocated memory inside this function. It is
 * the caller's responsibility to free this memory.
 */
gsl_vector* corr_norm_rvars(const int NUM_ASSETS, gsl_rng *rng,
        gsl_matrix *cholesky);

/* Calculates the return of a portfolio for one month given distribution
 * information for two risky assets and two correlated random variables
 *
 * \param assets an array of risky_asset structs that contains the statistics
 *  for each asset
 * \param NUM_ASSETS the number of assets in the assets array
 * \param rans a vector of correlated random variables
 * \return the percentage return for the portfolio over the course of a month
 */
double one_month_portfolio_return(risky_asset assets[],
        const int NUM_ASSETS, gsl_vector *rans);

/* Calculates the return of a single risky asset given distribution
 * information for that asset and a random variable for that month
 *
 * \param mean the annualized mean of the asset's returns
 * \param sigma the annualized standard deviation of the asset's returns
 * \param rand_var a random variable that represents the magnitude at which
 *  the asset's price increased or decreased over a one-month period
 * \return the percentage return on the single asset over a one-month period
 */
double one_month_asset_return(double mean, double sigma, double rand_var);
#endif
