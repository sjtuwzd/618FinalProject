//
// Created by zheng on 4/24/2020.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include <curl/curl.h>
#include <math.h>

char **read_ticker_file(char *filename, size_t *NUM_STOCKS) {

    FILE *ticker_file = fopen(filename, "r");
    char **tickers;

    const size_t  MAX_TICKER_SIZE = 6;
    if(ticker_file) {
        char line[MAX_TICKER_SIZE + 1];

        // GET THE FIRST LINE, WHICH IS THE TOTAL NUMBER
        fgets(line, sizeof(line), ticker_file);
        // sscanf return the number of parameters successfully filled
        int scanned = sscanf(line, "%zu", NUM_STOCKS);
        tickers = malloc(*NUM_STOCKS * sizeof(char *));
        for(int i = 0; i < *NUM_STOCKS && scanned != EOF; i++) {
            //get the whole line
            fgets(line, sizeof(line), ticker_file);
            char *temp = malloc((strlen(line) + 1) * sizeof(char));
            scanned = sscanf(line, "%5s", temp);
            tickers[i] = temp;
        }
    } else{
        printf("ERROR: Ticker file could not be opened");
        exit(1);
    }
    fclose(ticker_file);
    return tickers;
}



double *read_weight_file(char *filename, const size_t NUM_STOCKS){
    FILE *file = fopen(filename, "r");
    double *weights = malloc(NUM_STOCKS * sizeof(double ));
    if(file) {
        int scanned = 0;
        double sum = 0;
        for(int i = 0; i < NUM_STOCKS && scanned != EOF; i++) {
            // long double
            scanned = fscanf(file, "%lg", &weights[i]);
        }
        for (int j = 0; j < NUM_STOCKS ; ++j) {
            sum += weights[j];
        }
    } else {
        printf("ERROR: Weights file could not be opened\n");
        exit(1);
    }
    return weights;
}

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

    FILE *data_file = fopen(filename, "r");

    /* Monthly data for the past five years (~60 months) for each stock will
     * most likely be the number of data points used in each stock
     * price CSV file from Yahoo */
    size_t ret_cap = 72;
    size_t num_rets = 0;
    double *ret_data = malloc(ret_cap * sizeof(double));

    if (data_file) {
        /* Initialize to negative values to be able to tell whether the first
         * price has been read */
        double prev_price = -1, curr_price = -1;

        /* Throw away the first line of input
         * It just contains header labels */
        int scanned = fscanf(data_file, "%*s");

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
                                 "%*[^,],%*g,%*g,%*g,%*g,%*d,%lg", &curr_price);
            } else {
                scanned = fscanf(data_file,
                                 "%*[^,],%*g,%*g,%*g,%*g,%*d,%lg", &prev_price);
                if (scanned != EOF) {
                    /* Formula for monthly return =
                     * ln(this month's price/last month's price) */
                    ret_data[num_rets] = log(curr_price / prev_price);
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
    return ret_data;
}
