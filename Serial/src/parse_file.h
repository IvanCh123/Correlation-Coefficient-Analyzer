#ifndef PARSE_FILE_H
#define PARSE_FILE_H

#include <stdbool.h>

typedef struct
{		
	int row_count;
	int column_count;
	double** values;
	char** names;
	char** gens;
} data_t;


/**
    * @brief Fill in a matrix with the values stored in the specified CSV file.
    * @param input Name of the CSV file which contains the data set to be summarized. 
    * @param data  An struct containing the matrix to fill and it's dimensions.
    * @param transpose True if the given data set is transposed. 
    */
    
void load_file(char *input_file, data_t* data, bool transpose);


/**
    * @brief Free the memory required to store the data set given by the user.
    * @param data An struct containing the data set and it's dimensions.
    * @param transpose True if the given data set is transposed. 
    */
void parse_destroy(data_t* data, bool transpose);


#endif // CORRELATION_COEFFICIENT_SUMMARIZER_H
