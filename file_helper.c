//
// Created by zheng on 4/24/2020.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
