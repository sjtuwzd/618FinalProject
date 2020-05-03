#include "driver_lib.h"
#include <gsl/gsl_linalg.h>
#include <ctype.h>
#include <string.h>
#include <curl/curl.h>
#include <math.h>
#include <gsl/gsl_statistics_double.h>

#ifndef GSL_MAT_H
#define GSL_MAT_H
#include <gsl/gsl_matrix.h>
#endif

char *str_upr(char *upr, char *src) {
    /* Iterate through the source string until null character is reached */
    int final = 0;
    for (int i = 0; i < strlen(src); i++) {
        upr[i] = toupper(src[i]);
        final = i;
    }
    /* Terminate upr with the null character */
    upr[final + 1] = 0;
    return upr;
}

char *get_stock_file(char *ticker, struct tm end_date, int num_years) {

    char *tick_upr = malloc((strlen(ticker) + 1) * sizeof(char));
    /* Yahoo Finance URL uses upper case ticker as a parameter */
    str_upr(tick_upr, ticker);

    /* Build URL string for yahoo finance - URL usually contains ~ 100 chars,
     * so a 200 char buffer will suffice */
    size_t url_chars = 200;
    char url[url_chars];
    snprintf(url, url_chars, "http://real-chart.finance.yahoo.com/table.csv?"
            "s=%s&a=%02d&b=%02d&c=%04d&d=%02d&e=%02d&f=%04d&g=m&ignore=.csv",
            tick_upr, end_date.tm_mon, end_date.tm_mday,
            end_date.tm_year + 1900 - num_years, end_date.tm_mon,
            end_date.tm_mday, end_date.tm_year + 1900);

    /* Create the file name string, simply the stock ticker with an extension */
    char *ext = "csv";
    char *subdir = "data/prices";

    /* Allocate memory for the name to save the file under, create the filename
     * string, and open it for writing.
     * Add 3 extra characters for the '/', '.' and null characters */
    size_t num_chars = strlen(subdir) + strlen(tick_upr) + strlen(ext) + 3;
    char *filename = malloc(num_chars * sizeof(char));
    snprintf(filename, num_chars, "%s/%s.%s", subdir, tick_upr, ext);
    free(tick_upr);
    FILE *file = fopen(filename, "w");

    /* libcurl Calls */

    /* Initialize libcurl library objects */
    curl_global_init(CURL_GLOBAL_ALL);
    CURL *handle = curl_easy_init();
    /* Set the URL to request and file to be written to */
    curl_easy_setopt(handle, CURLOPT_URL, url);
    curl_easy_setopt(handle, CURLOPT_WRITEDATA, file);
    /* Perform curl action - request file from URL and save the file */
    CURLcode res = curl_easy_perform(handle);
    if (res != CURLE_OK)
        printf("LIBCURL ERROR: %s\n", curl_easy_strerror(res));
    /* Free resources taken up by libcurl library objects */
    curl_global_cleanup();

    fclose(file);

    return filename;
}

double *read_price_file(char *filename, size_t *data_size) {

    printf("Entering read_price_file of %s\n", filename);

    FILE *data_file = fopen(filename, "r");

    /* Monthly data for the past five years (~60 months) for each stock will
     * most likely be the number of data points used in each stock
     * price CSV file from Yahoo */
    size_t ret_cap = 72;
    size_t num_rets = 0;
    double *ret_data = malloc(ret_cap * sizeof(double));

    if (data_file) {
        printf("%s file exist, if passed!\n", filename);
        /* Initialize to negative values to be able to tell whether the first
         * price has been read */
        double prev_price = -1, curr_price = -1;

        /* Throw away the first line of input
         * It just contains header labels */
        int scanned = fscanf(data_file, "%*s");
        printf("%d", scanned);
        while (scanned != EOF) {
            /* Make sure that array is large enough to store the next
             * data point */
            if (num_rets >= ret_cap) {
                /* Double the size of the array more space is needed */
                ret_cap *= 2;
                double *tmp = realloc(ret_data, ret_cap * sizeof(double));
                if (!tmp) {
                    printf("ERROR: Array could not be expanded\n");
                    exit(1);
                }
                ret_data = tmp;
            }

            /* If this is the first line to be read, 
             * there is no current price, so read in first price */
            if (curr_price < 0) {
                scanned = fscanf(data_file,
                        "%*[^,],%*g,%*g,%*g,%*g,%lg,%*d", &curr_price);
            } else {
                scanned = fscanf(data_file,
                        "%*[^,],%*g,%*g,%*g,%*g,%lg,%*d", &prev_price);
                if (scanned != EOF) {
//                    printf("current price is %f, previous price is%f\n", curr_price, prev_price);
                    /* Formula for monthly return =
                     * ln(this month's price/last month's price) */
                    ret_data[num_rets] = log(curr_price / prev_price);
                    printf("%s: %f\n", filename, ret_data[num_rets]);
                    num_rets++;

                    /* Parallelism inhibitor: loop-carried dependency */
                    curr_price = prev_price;
                }
            }
        }
    } else {
        printf("ERROR: Could not open file\n");
        exit(1);
    }

    /* Set the value for the size of the data to be accessible by the caller */
    *data_size = num_rets;
    fclose(data_file);
    return ret_data;
}

char **read_ticker_file(char *filename, size_t *NUM_STOCKS) {
    FILE *ticker_file = fopen(filename, "r");
    char **tickers;
    
    /* Assuming the maximum size for a stock ticker is 6 characters.
     * Adjust this if need be. Keeping it low ensures that not too much
     * memory is wasted. */
    const size_t MAX_TICK_SIZE = 6;
    if (ticker_file) {
        /* Each line only contains one ticker in the tickers file */
        char line[MAX_TICK_SIZE + 1];

        /* Get the first line, which is just the number of stocks to be read */
        fgets(line, sizeof line, ticker_file);
        int scanned = sscanf(line, "%lu", NUM_STOCKS);

        tickers = malloc(*NUM_STOCKS * sizeof(char *));
        /* Read the tickers and populate the tickers array */
        for (int i = 0; i < *NUM_STOCKS && scanned != EOF; i++) {
            /* Get the whole line - safe way to prevent buffer overflow */
            fgets(line, sizeof line, ticker_file);
            /* Allocate memory to the current ticker string based on the size
             * of the line, then add it to the array */
            char *temp = malloc((strlen(line) + 1) * sizeof(char));
            scanned = sscanf(line, "%5s", temp);
            tickers[i] = temp;
        }
    } else {
        printf("ERROR: Ticker file could not be opened\n");
        exit(1);
    }
    fclose(ticker_file);
    
    return tickers;
}

double *read_weight_file(char *filename, const size_t NUM_STOCKS) {
    FILE *file = fopen(filename, "r");
    double *weights = malloc(NUM_STOCKS * sizeof(double));
    if (file) {
        int scanned = 0;
        double sum = 0;
        /* Get the current weight and add to array */
        for (int i = 0; i < NUM_STOCKS && scanned != EOF; i++)
            scanned = fscanf(file, "%lg", &weights[i]);
        /* Total up all weights in the array and check if it is equal to 100% */
        for (int i = 0; i < NUM_STOCKS; i++)
            sum += weights[i];
    } else {
        printf("ERROR: Weights file could not be opened\n");
        exit(1);
    }
    return weights;
}

gsl_matrix *calculate_varcovar(ret_data *dataset, size_t NUM_STOCKS) {

    /* Create the varcorvar matrix */
    gsl_matrix *varcovar = gsl_matrix_alloc(NUM_STOCKS, NUM_STOCKS);

    /* Fill the upper triangle (and diagonal) of the varcovar matrix 
     * by calculating the covariance */
    for (int i = 0; i < NUM_STOCKS; i++) {
        for (int j = i; j < NUM_STOCKS; j++) {
            /* Assumes that datasets will have the same size - a reasonable
             * assumption since all stocks will have return data going back
             * to and from the same points in time */
            double covar = gsl_stats_covariance(dataset[i].data, 1,
                    dataset[j].data, 1, dataset[i].size);
            printf("Current i is %d, current j is %d, covar is %f\n",i, j, covar);
            gsl_matrix_set(varcovar, i, j, covar);
        }
    }

    /* Fill lower triangle of the varcovar matrix */
    for (int i = 1; i < NUM_STOCKS; i++) {
        for (int j = 0; j < i; j++) {
            gsl_matrix_set(varcovar, i, j,
                    gsl_matrix_get(varcovar, j, i));
        }
    }
    return varcovar;
}

gsl_matrix* varcovar_from_file(const char *filename, int *NUM_ASSETS) {
    gsl_matrix *varcovar;
    FILE *vc_file = fopen(filename, "r");
    if (vc_file) {

        /* Get the number of assets in the portfolio as first line in file */
        fscanf(vc_file, "%d", NUM_ASSETS);

        /* Create matrix with size NUM_ASSETS by NUM_ASSETS */
        varcovar = gsl_matrix_alloc(*NUM_ASSETS, *NUM_ASSETS);

        /* Read all numbers from file and place into matrix */
        for (int i = 0; i < *NUM_ASSETS; i++) {
            for (int j = 0; j < *NUM_ASSETS; j++) {
                double curr; 
                fscanf(vc_file, "%lf,", &curr);
                gsl_matrix_set(varcovar, i, j, curr);
            } 
        }
    } else {
        printf("ERROR: no variance-covariance matrix file provided\n");
        exit(1);
    }
    fclose(vc_file);
    return varcovar;
}

risky_asset* assets_from_file(const char *filename,
        const int NUM_ASSETS) {

    risky_asset *assets = (risky_asset *)malloc(NUM_ASSETS *
            sizeof(risky_asset));

    FILE *asset_file = fopen(filename, "r");
    if (asset_file) {
        for (int i = 0; i < NUM_ASSETS; i++) {
            fscanf(asset_file, "%lg,%lg,%lg", &assets[i].mean, &assets[i].sigma,
                    &assets[i].port_weight);
        }
    } else {
        printf("ERROR: valid asset file not provided\n");
        exit(1);
    }
    fclose(asset_file);

    /* Ensure the weight of the whole portfolio is equal to 100% */
    double total_weight = 0;
    for (int i = 0; i < NUM_ASSETS; i++) {
        total_weight += assets[i].port_weight;
    }
    if (total_weight != 1.0) {
        printf("ERROR: weights of assets do not total 100%%\n");
        exit(1);
    }

    return assets;
}


