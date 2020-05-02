#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "file_helper.h"
#include <string.h>
#include "portfolio_helper.h"
#include <time.h>
//#include <curl.h>
#include <sys/time.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_matrix.h>

int main(int argc, char **argv) {
    printf("Hello, World!\n");
    unsigned long seed = 0;

    if(argc < 4 || argc > 5) {
        printf("ERROR: ARGUMENTS NUMBER INVALID");
        exit(1);
    }
    else if(argc == 5) {
        seed = strtoul(argv[4], NULL, 10);
    }

    // get the input arguments
    const int NUM_MONTHS = atoi(argv[1]);
    const int NUM_RUNS = atoi(argv[2]);
    const char *PORT_SUFF = argv[3];
    size_t NUM_ASSETS = 0;

    char *subdir = "data";
    char *tickers_name = "tickers";
    char *weights_name = "weights";
    char *extension = "csv";
    char ticker_filename[strlen(subdir) + strlen(tickers_name)
                         + strlen(PORT_SUFF) + strlen(extension) + 3];
    char weights_filename[strlen(subdir) + strlen(weights_name)
                          + strlen(PORT_SUFF) + strlen(extension) + 3];

    // write formatted output to sized buffer
    snprintf(ticker_filename, sizeof(ticker_filename), "%s/%s%s.%s",
             subdir, tickers_name, PORT_SUFF, extension);
    snprintf(weights_filename, sizeof(weights_filename), "%s/%s%s.%s",
             subdir, weights_name, PORT_SUFF, extension);


    char **ticks = read_ticker_file(ticker_filename, &NUM_ASSETS);
    double *weights = read_weight_file(weights_filename, NUM_ASSETS);

    /* Test if the files are correctly read*/
    for (int i = 0; i < NUM_ASSETS; ++i) {
        printf("%s\n", ticks[i]);
        printf("%lg\n", weights[i]);
    }


    ret_data *dataset = malloc(NUM_ASSETS * sizeof(ret_data));
    risky_asset *assets = malloc(NUM_ASSETS * sizeof(risky_asset));

    time_t t = time(NULL);
    struct tm curr_time = *localtime(&t);

    /* Set up data structs for measuring the time taken to both retrieve
  * stock data and run the simulations */
    struct timeval retrieval_begin, retrieval_end, sim_begin, sim_end;
    double retrieval_time, sim_time;
    gettimeofday(&retrieval_begin, NULL);

//    for(int i = 0; i < NUM_ASSETS; i++) {
//        char*
//    }

    for(int i = 0; i < NUM_ASSETS; i++) {
//        char *price_filename = get_stock_file(ticks[i], curr_time, 6);
//        dataset[i].data = read_price_file(price_filename, &dataset[i].size);

        /* Read the price data from locally-stored files */

        char *price_subdir = "data/prices";
        char *price_extension = "csv";
        size_t price_fname_chars = strlen(ticks[i]) + strlen(price_subdir) +
                        strlen(price_extension) + 3;
        char price_filename[price_fname_chars];
        snprintf(price_filename, price_fname_chars, "%s/%s.%s",
                price_subdir, ticks[i], price_extension);
        dataset[i].data = read_price_file(price_filename, &dataset[i].size);


        assets[i].ticker = malloc((strlen(ticks[i]) + 1) * sizeof(char));
        strcpy(assets[i].ticker, ticks[i]);

        assets[i].mean = gsl_stats_mean(dataset[i].data,1,dataset[i].size) * 12;
        assets[i].sigma = gsl_stats_sd(dataset[i].data,1,dataset[i].size) * sqrt(12);
        assets[i].port_weight = weights[i];
    }

    for (int i = 0; i < NUM_ASSETS; i++) {
        free(ticks[i]);
    }

    free(ticks);
    free(weights);

    /* Calculate the variance-covariance based on the price datasets of
  * all assets */
    gsl_matrix *varcovar = calculate_varcovar(dataset, NUM_ASSETS);
    free(dataset);
    perform_cholesky(varcovar, NUM_ASSETS);
    gsl_matrix *cholesky = varcovar;

    printf("Finished crunching numbers and getting data for stock returns\n"
           "Beginning simulations\n");

    /* This file will contain the final porfolio returns for all runs */
    FILE *results_file = fopen("data/results.txt","w");
    double total_rets[NUM_RUNS];

    /* Start the Monte Carlo simulation timer */
    gettimeofday(&sim_begin, NULL);

    if(results_file) {
        /* Each run of this loop is one simulation of portfolio's return */
        for (int run = 0; run < NUM_RUNS; run++) {
            double total_return = 0;

            gsl_rng *rng;
            /* Seed is shared, so lock it down when it is being used and
             * incremented so that no two RNGs are initialized with the
             * same seed */
            {
                rng = initialize_rng_with_seed(seed);
                /* When seed hits max unsigned long it will wrap to 0 */
                seed++;
            }

            /* Calculate the monthly portfolio return
             * and add to yearly return accumulator */
            for (int month = 1; month <= NUM_MONTHS; month++) {
                gsl_vector *rans = corr_norm_rvars(NUM_ASSETS, rng, cholesky);
                double month_ret = one_month_portfolio_return(assets,
                                                              NUM_ASSETS, rans);
                gsl_vector_free(rans);
                total_return += month_ret;
            }
            gsl_rng_free(rng);

            total_rets[run] = total_return;
        }
        gettimeofday(&sim_end, NULL);

        /* Calculate time elapsed for data retrieval and simulations */
        retrieval_time = retrieval_end.tv_sec - retrieval_begin.tv_sec;
        retrieval_time += ((double)(retrieval_end.tv_usec -
                                    retrieval_begin.tv_usec)) / 1000000.0;

        sim_time = sim_end.tv_sec - sim_begin.tv_sec;
        sim_time += ((double)(sim_end.tv_usec - sim_begin.tv_usec)) / 1000000.0;
        printf("Data retrieval time: %lg\n"
               "Simulation time: %lg\n", retrieval_time, sim_time);

        for (int i = 0; i < NUM_RUNS; i++) {
            fprintf(results_file, "%lg\n", total_rets[i] * 100);
        }
        fclose(results_file);

        double res_mean = gsl_stats_mean(total_rets,1,NUM_RUNS) * 100;
        double res_std = gsl_stats_sd(total_rets,1,NUM_RUNS)* 100;

        printf("Mean of annual returns: %lg%%\n"
               "Standard dev of annual returns: %lg%%\n",
               res_mean, res_std);
    } else {
        printf("NO RESULTS FILE\n");
    }
    gsl_matrix_free(cholesky);
    return 0;
}



