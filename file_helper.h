//
// Created by zheng on 4/24/2020.
//

#ifndef INC_618FINALPROJECT_FILE_HELPER_H
#define INC_618FINALPROJECT_FILE_HELPER_H

#endif //INC_618FINALPROJECT_FILE_HELPER_H


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