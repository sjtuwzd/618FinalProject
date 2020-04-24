#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "file_helper.h"
#include <string.h>
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

    return 0;
}
