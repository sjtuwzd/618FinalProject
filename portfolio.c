#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>

#include <omp.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics_double.h>

#include "portfolio_lib.h"
#include "driver_lib.h"

#ifndef MPI
#define MPI 0
#endif

#if MPI
#include <mpi.h>
void gather_perf(int nzone, double *total_rets, int NUM_RUNS);
void send_perf(double *total_rets, int run_start, int run_end);
void gather_node_state(int nzone, char* add_ticker_0, char* decrease_ticker_0, double best_result_0, double best_dev_0);
void send_node_state(char* add_ticker, char* decrease_ticker, double best_result, double best_dev);
#endif

int main(int argc, char **argv) {

    unsigned long seed = 0;
    if (argc < 5 || argc > 6) {
        printf("ERROR: You have not provided exactly five or "
               "exactly six arguments\n");
        exit(1);
        /* Set the seed if it has been provided */
    } else if (argc == 6) {
        seed = strtoul(argv[5], NULL, 10);
    }

    const int NUM_MONTHS = atoi(argv[1]);
    const int NUM_RUNS = atoi(argv[2]);
    const char *PORT_SUFF = argv[3];
    char *METHOD = argv[4];
    printf("%s", METHOD);
    size_t NUM_ASSETS = 0;

    /* Build the names of the files containing the assets' tickers and
     * their weights in the portfolio */
    char *subdir = "data";
    char *tickers_name = "tickers";
    char *weights_name = "weights";
    char *extension = "csv";
    char ticker_filename[strlen(subdir) + strlen(tickers_name)
                         + strlen(PORT_SUFF) + strlen(extension) + 3];
    char weights_filename[strlen(subdir) + strlen(weights_name)
                          + strlen(PORT_SUFF) + strlen(extension) + 3];

    /* Takes format: 'data/tickersX.csv' and 'data/weightsX.csv' where X is
     * the suffix that identifies the specific portfolio */
    snprintf(ticker_filename, sizeof(ticker_filename), "%s/%s%s.%s",
             subdir, tickers_name, PORT_SUFF, extension);
    snprintf(weights_filename, sizeof(weights_filename), "%s/%s%s.%s",
             subdir, weights_name, PORT_SUFF, extension);

    /* Read file into an array of ticker strings */
    char **ticks = read_ticker_file(ticker_filename, &NUM_ASSETS);
    /* Read file into an array of weight percentages (in decimal form) */
    double *weights = read_weight_file(weights_filename, NUM_ASSETS);

    ret_data *dataset = malloc(NUM_ASSETS * sizeof(ret_data));
    risky_asset *assets = malloc(NUM_ASSETS * sizeof(risky_asset));

    time_t t = time(NULL);
    struct tm curr_time = *localtime(&t);

    /* Set up data structs for measuring the time taken to both retrieve
     * stock data and run the simulations */
    struct timeval retrieval_begin, retrieval_end, sim_begin, sim_end, total_sim_begin, total_sim_end;
    double retrieval_time, sim_time, total_sim_time;

    /* Start data retrieval timer */
    gettimeofday(&retrieval_begin, NULL);
    // #pragma omp parallel for num_threads(NUM_ASSETS)
    for (int i = 0; i < NUM_ASSETS; i++) {
        printf("start reading price file!\n");
        /* Fetch price data from Yahoo! Finance */
        /* char *price_filename = get_stock_file(ticks[i], curr_time, 6);
         dataset[i].data = read_price_file(price_filename, &dataset[i].size);*/
        /* Read the price data from locally-stored files */

        char *price_subdir = "data/prices";
        char *price_extension = "csv";
        size_t price_fname_chars = strlen(ticks[i]) + strlen(price_subdir) +
                                   strlen(price_extension) + 3;
        char price_filename[price_fname_chars];
        snprintf(price_filename, price_fname_chars, "%s/%s.%s",
                 price_subdir, ticks[i], price_extension);

        printf("%s\n", price_filename);
        dataset[i].data = read_price_file(price_filename, &dataset[i].size);


        assets[i].ticker = malloc((strlen(ticks[i]) + 1) * sizeof(char));

        /* Actually copy the string instead of setting pointers equal so that
         * the ticks array can be freed */
        strcpy(assets[i].ticker, ticks[i]);

        /* Set the rest of the current asset's attributes */
        // assets[i].mean = gsl_stats_mean(dataset[i].data, 1, dataset[i].size) * 12;
        assets[i].mean = gsl_stats_mean(dataset[i].data, 1, dataset[i].size) * 252;
        // assets[i].sigma = gsl_stats_sd(dataset[i].data, 1, dataset[i].size) * sqrt(12);
        assets[i].sigma = gsl_stats_sd(dataset[i].data, 1, dataset[i].size) * sqrt(252);
        assets[i].port_weight = weights[i];
        printf("end reading price file!\n");
    }
    /* End data retrieval timer */
    gettimeofday(&retrieval_end, NULL);

    for (int i = 0; i < NUM_ASSETS; i++)
        free(ticks[i]);
    free(ticks);
    free(weights);

    printf("End freeing memory!\n");
    /* Calculate the variance-covariance based on the price datasets of
     * all assets */
    gsl_matrix *varcovar = calculate_varcovar(dataset, NUM_ASSETS);

    // for (int i = 0; i < 2; i ++) {
    //     for (int j = 0; j < 2; j++) {
    //         printf("%f\n", gsl_matrix_get(varcovar, i, j));
    //     }
    // }

    printf("End calculating varcovar!\n");
    free(dataset);
    printf("Entering perform_cholesky!\n");
    perform_cholesky(varcovar, NUM_ASSETS);
    printf("Leaving perform_cholesky!\n");
    gsl_matrix *cholesky = varcovar;

    printf("Finished crunching numbers and getting data for stock returns\n"
           "Beginning simulations\n");

    double change_amount = 0.05;
    // variables used to check the best strategy
    double best_result = -100.0;
    int best_i = 0;
    int best_j = 0;
    double best_dev = 0.0;

    bool best_mode = false;

    printf("\nStart MPI.\n");

    int process_count = 1;
    int this_zone = 0;
    int nzone = 0;
#if MPI
    if (strcmp(METHOD, "best") == 0) {
        best_mode = true;
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &process_count);
        MPI_Comm_rank(MPI_COMM_WORLD, &this_zone);
    }
#endif
    int i_start, i_end;
    int j_end = NUM_ASSETS;
    bool mpi_master = this_zone == 0;

    // if is for best strategy
    if (strcmp(METHOD, "best") == 0) {
        nzone = process_count;
        printf("Rank %d out of %d processes.\n", this_zone, nzone);

        i_start = this_zone * (ceil((double)NUM_ASSETS / (double)nzone));
        i_end = (this_zone + 1) * (ceil((double)NUM_ASSETS / (double)nzone));
        i_end = this_zone == nzone - 1 ? NUM_ASSETS : i_end;

        seed = this_zone * NUM_ASSETS * NUM_ASSETS * NUM_RUNS;
        // for (int i = 0; i < NUM_ASSETS; i++) {
    }
    // just for predicting the current performance
    else {
        printf("Hello!\n");
        i_start = 0;
        i_end = 1;
        j_end = 1;
    }

    // start timing
    if (mpi_master) {
        gettimeofday(&total_sim_begin, NULL);
    }
    printf("J end is: %d\n", j_end);
    printf("Process %d, start i is %d, end i is %d\n", this_zone, i_start, i_end);
    for (int i = i_start; i < i_end; i++) {
        for (int j = 0; j < j_end; j++) {

            // not the same stock or the first run
            if (i == 0 || i != j) {
                double original_i_weight = assets[i].port_weight;
                double original_j_weight = assets[j].port_weight;


                if (assets[i].port_weight >= change_amount) {
                    assets[i].port_weight -= change_amount;
                    assets[j].port_weight += change_amount;
                } else {
                    assets[i].port_weight = 0.0;
                    assets[j].port_weight += change_amount;
                }
                char *result_subdir = "data";
                char *result_extension = "txt";
                char *name = "result";
                size_t result_fname_chars = strlen(result_subdir) + strlen(assets[i].ticker) + strlen(assets[j].ticker) +
                                            strlen(result_extension) + 13 + strlen(name);
                char result_filename[result_fname_chars];

                snprintf(result_filename, result_fname_chars, "%s/%s_add%s_minus%s.%s",
                         result_subdir, name, assets[j].ticker, assets[i].ticker, result_extension);

                // printf("Adding weight of %s, decreasing weight of %s\n", assets[j].ticker, assets[i].ticker);
                /* This file will contain the final porfolio returns for all runs */
                // FILE *results_file = fopen("data/results.txt", "w");
                FILE *results_file = fopen(result_filename,  "w");
                double total_rets[NUM_RUNS];

                /* Start the Monte Carlo simulation timer */
                gettimeofday(&sim_begin, NULL);

                if (results_file) {
                    int run_start = 0, run_end = NUM_RUNS;
                    /* Each run of this loop is one simulation of portfolio's return */
                    // #pragma omp parallel for
#if MPI
                    if (!best_mode) {
                        MPI_Init(NULL, NULL);
                        MPI_Comm_size(MPI_COMM_WORLD, &process_count);
                        MPI_Comm_rank(MPI_COMM_WORLD, &this_zone);

                        nzone = process_count;
                        mpi_master = this_zone == 0;
                        printf("Performance Simulation: rank %d out of %d processes.\n", this_zone, nzone);

                        run_start = this_zone * (ceil((double)NUM_RUNS / (double)nzone));
                        run_end = (this_zone + 1) * (ceil((double)NUM_RUNS / (double)nzone));
                        run_end = this_zone == nzone - 1 ? NUM_RUNS : run_end;
                        printf("run_start = %d, run_end = %d\n", run_start, run_end);
                    }
#endif

                    for (int run = run_start; run < run_end; run++) {
                        double total_return = 0;

                        gsl_rng *rng;
                        /* Seed is shared, so lock it down when it is being used and
                         * incremented so that no two RNGs are initialized with the
                         * same seed */
                        // #pragma omp critical
                        // {
                        rng = initialize_rng_with_seed(seed);
                        /* When seed hits max unsigned long it will wrap to 0 */
                        seed++;
                        // }

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

#if MPI
                    if (!best_mode) {
                        if (mpi_master)
                            gather_perf(nzone, total_rets, NUM_RUNS);
                        else
                            send_perf(total_rets, run_start, run_end);

                        MPI_Finalize();
                    }
#endif

                    /* Calculate time elapsed for data retrieval and simulations */
                    retrieval_time = retrieval_end.tv_sec - retrieval_begin.tv_sec;
                    retrieval_time += ((double) (retrieval_end.tv_usec -
                                                 retrieval_begin.tv_usec)) / 1000000.0;

                    sim_time = sim_end.tv_sec - sim_begin.tv_sec;
                    sim_time += ((double) (sim_end.tv_usec - sim_begin.tv_usec)) / 1000000.0;
                    printf("Data retrieval time: %lg\n"
                           "Simulation time: %lg\n", retrieval_time, sim_time);

                    for (int run = 0; run < NUM_RUNS; run++) {
                        fprintf(results_file, "%lg\n", total_rets[run] * 100);
                    }
                    fclose(results_file);

                    double res_mean = gsl_stats_mean(total_rets, 1, NUM_RUNS) * 100;
                    double res_std = gsl_stats_sd(total_rets, 1, NUM_RUNS) * 100;

                    // the first run, store into the variable
                    if (i == 0 && j == 0) {
                        best_result = res_mean;
                        best_dev = res_std;
                        best_i = i;
                        best_j = j;
                    } else {
                        if (res_mean > best_result) {
                            best_result = res_mean;
                            best_dev = res_std;
                            best_i = i;
                            best_j = j;
                        }
                    }
                    if (best_mode || (!best_mode && mpi_master))
                        printf("Process %d, I is %d, j is %d, Adding weight of %s, decreasing weight of %s, Mean of percentage returns: %lg%%.Standard dev of percentage returns: %lg%%\n", this_zone, i, j,
                               assets[j].ticker, assets[i].ticker, res_mean, res_std);
                } else {
                    printf("NO RESULTS FILE\n");
                }

                assets[i].port_weight = original_i_weight;
                assets[j].port_weight = original_j_weight;
            }

        }
    }

    // if (mpi_master) {
    //     printf("Process %d: The best strategy is to add %s while decreasing %s\n", this_zone, assets[best_j].ticker, assets[best_i].ticker);
    //     printf("Process %d: The Resulted Mean of percentage returns: %lg%%\n"
    //            "The Resulted Standard dev of percentage returns: %lg%%\n",
    //            this_zone, best_result, best_dev);
    // }

#if !MPI
    printf("Process %d: The best strategy is to add %s while decreasing %s\n", this_zone, assets[best_j].ticker, assets[best_i].ticker);
    printf("Process %d: The Resulted Mean of percentage returns: %lg%%\n"
           "The Resulted Standard dev of percentage returns: %lg%%\n",
           this_zone, best_result, best_dev);
#endif
    gsl_matrix_free(cholesky);

#if MPI
    if (strcmp(METHOD, "best") == 0) {
        if (mpi_master)
            gather_node_state(nzone, assets[best_j].ticker, assets[best_i].ticker, best_result, best_dev);
        else
            send_node_state(assets[best_j].ticker, assets[best_i].ticker, best_result, best_dev);

        MPI_Finalize();
    }
#endif

    if (mpi_master) {
        gettimeofday(&total_sim_end, NULL);
        total_sim_time = total_sim_end.tv_sec - total_sim_begin.tv_sec;
        total_sim_time += ((double) (total_sim_end.tv_usec - total_sim_begin.tv_usec)) / 1000000.0;
        printf("Total Simulation time: %lg\n", total_sim_time);
    }
    exit(0);
}

#if MPI
/* Called by process 0 to collect performance states from all other processes */
void gather_perf(int nzone, double *total_rets, int NUM_RUNS) {
    int i, j, run_start, run_end, run_num;

    for (i = 1; i < nzone; i++) {
        run_start = i * (ceil((double)NUM_RUNS / (double)nzone));
        run_end = (i + 1) * (ceil((double)NUM_RUNS / (double)nzone));
        run_end = i == nzone - 1 ? NUM_RUNS : run_end;
        run_num = run_end - run_start;
        printf("Receive %d total return valus from process %d\n", run_num, i);

        double *data_buf = malloc(run_num * sizeof(double));

        MPI_Recv(data_buf, run_num,
                 MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        for (j = 0; j < run_num; j++) {
            total_rets[run_start + j] = data_buf[j];
        }

        free(data_buf);
    }
}

/* Called by other processes to send their node states to process 0 */
void send_perf(double *total_rets, int run_start, int run_end) {
    int run_num = run_end - run_start;

    MPI_Send(total_rets + run_start, run_num,
             MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

/* Called by process 0 to collect node states from all other processes */
void gather_node_state(int nzone, char* add_ticker_0, char* decrease_ticker_0, double best_result_0, double best_dev_0) {
    int i;
    /* MPI buffer */
    // char add_ticker[5];
    // char decrease_ticker[5];
    // double data_buf[2];

    int curr_best_zone = 0;
    char* curr_add_ticker = add_ticker_0;
    char* curr_decrease_ticker = decrease_ticker_0;
    double curr_best_result = best_result_0;
    double curr_best_dev = best_dev_0;

    for (i = 1; i < nzone; i++) {

        printf("At beginning decrease ticker is %s\n", curr_decrease_ticker);
        /* MPI buffer */
        // char add_ticker[5];
        // char decrease_ticker[5];
        char *add_ticker = malloc(5 * sizeof(char));
        char *decrease_ticker = malloc(5 * sizeof(char));
        double data_buf[2];

        MPI_Recv(add_ticker, 5,
                 MPI_CHAR, i, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        MPI_Recv(decrease_ticker, 5,
                 MPI_CHAR, i, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        MPI_Recv(data_buf, 2,
                 MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        printf("Received from process %d: The best strategy is to add %s and decrease %s.\n"
               "\tThe Resulted Mean of percentage returns: %lg%%\n"
               "\tThe Resulted Standard dev of percentage returns: %lg%%\n",
               i, add_ticker, decrease_ticker, data_buf[0], data_buf[1]);

        if (curr_best_result < data_buf[0]) {
            curr_best_zone = i;
            printf("Discovered %d is better! Exchanging!\n", i);
            curr_add_ticker = add_ticker;
            curr_decrease_ticker = decrease_ticker;
            curr_best_result = data_buf[0];
            curr_best_dev = data_buf[1];
        }
        printf("Current decrease ticker is %s\n", curr_decrease_ticker);

        free(add_ticker);
        free(decrease_ticker);
    }

    printf("\nThe best strategy is to add %s while decreasing %s, the best process is %d \n", curr_add_ticker, curr_decrease_ticker, curr_best_zone);
    printf("The Resulted Mean of percentage returns: %lg%%\n"
           "The Resulted Standard dev of percentage returns: %lg%%\n",
           curr_best_result, curr_best_dev);
}

/* Called by other processes to send their node states to process 0 */
void send_node_state(char* add_ticker, char* decrease_ticker, double best_result, double best_dev) {
    /* MPI buffer */
    double data_buf[2] = {best_result, best_dev};

    MPI_Send(add_ticker, 5,
             MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    MPI_Send(decrease_ticker, 5,
             MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    MPI_Send(data_buf, 2,
             MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}
#endif