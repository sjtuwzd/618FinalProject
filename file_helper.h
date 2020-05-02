//
// Created by zheng on 4/24/2020.
//

#ifndef INC_618FINALPROJECT_FILE_HELPER_H
#define INC_618FINALPROJECT_FILE_HELPER_H

#endif //INC_618FINALPROJECT_FILE_HELPER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include <curl/curl.h>
#include <math.h>
#include <gsl/gsl_matrix_double.h>
#include "portfolio_helper.h"

/**
 * Helper to read ticker file
 * @param filename
 * @param NUM_STOCKS
 * @return
 */
char **read_ticker_file(char *filename, size_t *NUM_STOCKS);


/**
 * Helper to read weight file
 * @param filename
 * @param NUM_STOCKS
 * @return
 */
double *read_weight_file(char *filename, const size_t NUM_STOCKS);



double *read_price_file(char *filename, size_t *data_size);


char *get_stock_file(char *ticker, struct tm end_date, int num_years);


gsl_matrix *calculate_varcovar(ret_data *dataset, size_t NUM_STOCKS);


gsl_matrix* varcovar_from_file(const char *filename, int *NUM_ASSETS);


risky_asset* assets_from_file(const char *filename,
                              const int NUM_ASSETS);