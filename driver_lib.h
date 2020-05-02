#ifndef GSL_MAT_H
#define GSL_MAT_H
#include <gsl/gsl_matrix.h>
#endif

#ifndef DRIVER_DEFS_H
#define DRIVER_DEFS_H

#include "portfolio_lib.h"
#include <time.h>

/* A struct to represent a set of return percentage data. It simply holds an
 * array of double data and the size of that data */
struct {
    double *data;
    size_t size;
} typedef ret_data;

/* Converts a string to upper case.
 *
 * \param upr is the buffer that will contain the upper case characters
 * \param src is the string whose characters will be converted to upper case
 *
 * Notes:
 * You must make sure that upr buffer is the same size as the src buffer before
 * calling this function. No memory is allocated inside this function.
 * */
char *str_upr(char *upr, char *src);

/* Retrieves a file of monthly prices for a given stock from Yahoo Finance
 * 
 * \param ticker is a string of the stock's ticker
 * \param start_time is a struct containing the ending date that the price data
 *  should be retrieved for
 * \param num_years is the number of years the price data will be retrieved for,
 *  ending with the end_date
 * \return the filename of the price file that was just retrieved
 *
 * Notes:
 * The memory for return value filename is allocated inside of this function.
 * The caller must free the memory taken by the returned value filename pointer.
 * There is no error checking for when the file requested is not returned by
 *  Yahoo.
 */
char *get_stock_file(char *ticker, struct tm end_date, int num_years);

/* Reads a file of monthly stock prices and calculates the percentage
 * return for each month
 *
 * \param filename is the name of the price file to be read
 * \param data_size is a pointer to the size of the return percentage dataset
 *  The value pointed to by data_size is modified in this function.
 * \return a pointer to an array of monthly return percentages
 *
 * Notes:
 * The memory for the returned array pointer is allocated in this function.
 * It is the caller's responsibility to free that memory.
 */
double *read_price_file(char *filename, size_t *data_size);

/* Reads a file containing the number of tickers to be read and the tickers of
 * stocks.
 *
 * \param filename is the name of the file to be read
 * \param NUM_STOCKS is a pointer to the number of stocks to be read. This
 *  function modifies the value pointed to.
 * \return an array of the ticker strings read
 *
 * Notes:
 * The memory for the returned array pointer and the ticker 
 * strings in that array is allocated in the function. The caller is
 * responsible for freeing all of that memory taken up both by the individual
 * strings and the array pointer itself.
 */
char **read_ticker_file(char *filename, size_t *NUM_STOCKS);

/* Read a file containing weights of stocks in a portfolio and create an
 * array of those weights
 *
 * \param filename is the name of the weight file to be read
 * \param NUM_STOCKS is the number of stocks whose weights will be read
 * \return a double array containing the weights of the stocks
 *
 * Notes:
 * The memory of the double array returned is allocated its memory inside
 * this function. The caller is responsbile for freeing all of that memory.
 */
double *read_weight_file(char *filename, const size_t NUM_STOCKS);

/* Calculates a variance-covariance matrix for a set of stocks whose returns
 * datasets are given
 *
 * \param dataset is an array of return datasets, and there is one dataset
 *  per stock in the portfolio
 * \param NUM_STOCKS is the number of stocks in the portfolio
 * \return a pointer to a variance-covariance matrix
 * 
 * Notes:
 * The memory of the gsl_matrix that is returned is allocated inside
 * this function. The caller is responsible for freeing all that memory.
 */
gsl_matrix *calculate_varcovar(ret_data *dataset, size_t NUM_STOCKS);

/* Initializes a variance covariance matrix based on the contents of a file
 *
 * \param filename is the name of the file with the variance covariance matrix
 *  data
 * \param NUM_ASSETS is a pointer to the number of stocks in the portfolio
 *  that the variance covariance matrix is for
 * \return a pointer to the variance covariance matrix
 * 
 * Notes:
 * This allocates memory for the matrix. It is the caller's responsibility to
 * free this memory. This function will also overwrite the value that NUM_ASSETS
 * points to
 */
gsl_matrix* varcovar_from_file(const char *filename, int *NUM_ASSETS);

/* Create an array of risky_assets from the contents of a file that specifies
 * the statistics of the distributions of those assets
 *
 * \param filename is the name of the file to be read
 * \param NUM_ASSETS is the number of assets to be read
 * \return an array of initialized assets
 *
 * Notes:
 * The returned array is allocated memory inside this function. It is the
 * caller's responsibility to free this memory
 */
risky_asset* assets_from_file(const char *filename,
        const int NUM_ASSETS);


#endif
