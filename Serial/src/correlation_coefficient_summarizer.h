#ifndef CORRELATION_COEFFICIENT_SUMMARIZER_H
#define CORRELATION_COEFFICIENT_SUMMARIZER_H

#include  "mathematical_operations.h"

struct data_set_info_t;
typedef struct data_set_info_t data_set_info_t;

/**
    * @brief Summarize the given data subset using as a discard mechanism the correlation or anticorrelation between its variables (cancer types).
    * @param data_set Complete data set to be reduced.
    * @param correlation_coeficients A matrix whose rows and columns associate two variables (cancer types) with their corresponding correlation coefficient.
    * @param variable_count Number of variables (cancer types).
    * @param subset_size Number of observations for each variable (gen types).
    * @param lower_bound Lower bound to determinate whether two variables are correlated-anticorrelated or not. 
    * @param upper_bound Upper bound to determinate whether two variables are correlated-anticorrelated or not. 
    * @param correlation_record An array whose cells correspond to each cancer type. It'll be one if the corresponding cancer type must be conserved, zero otherwise.
    * @param matches An array whose cells correspond to each cancer type. It'll be one if the corresponding cancer type matches the wildcard given by the user.
    */
void start_summarazing(double*** data_set, double*** correlation_coeficients,const size_t variable_count, const size_t subset_size, const double lower_bound, const double upper_bound, int** correlation_record, int* matches,const bool data_set_transposed);



#endif // CORRELATION_COEFFICIENT_SUMMARIZER_H
