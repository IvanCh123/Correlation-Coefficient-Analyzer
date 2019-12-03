#ifndef MATHEMATICAL_OPERATIONS_H
#define MATHEMATICAL_OPERATIONS_H

#include <stdbool.h>
#include <stddef.h>

/**
    * @brief Determine wheter two variables are correlated-anticorrelated or not.
    * @param pearson_correlation_coefficient Correlation coefficient between two variables.
    * @param lower_bound Lower bound to determinate whether two variables are correlated-anticorrelated or not. 
    * @param upper_bound Upper bound to determinate whether two variables are correlated-anticorrelated or not. 
    * @return True if given correlation coefficient is between the lower and upper bound.
    */
bool is_correlated(const double pearson_correlation_coefficient, const double lower_bound, const double upper_bound);

/**
    * @brief Calculates Pearson's correlation coefficient between two given variable subsets.
    * @param X_subset First variable observations. 
    * @param Y_subset Second variable observations.
    * @param subset_size Number of observations.
    * @return Pearson's correlation coefficient between two given variable subsets.
    */
double calculate_pearson_coeffient(double** X_subset, double** Y_subset, const size_t subset_syze);



#endif // MATHEMATICAL_OPERATIONS_H
