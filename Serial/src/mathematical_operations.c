#include <math.h>

#include "mathematical_operations.h"


/***
    * @brief Calculate the covariance between two given variable subsets.
    * @param X_subset First variable observations.
    * @param X_mean First variable mean. 
    * @param Y_subset Second variable observations.
    * @param Y_mean Second variable mean. 
    * @param subset_size Number of observations.
    * @return The covariance between two given variable subsets.
    */
double calculate_covariance(double** X_subset, const double X_mean, double** Y_subset, const double Y_mean, const size_t subset_syze);

/***
    * @brief Calculate the mean of a given variable subset.
    * @param subset Variable observations.
    * @param subset_size Number of observations.
    * @return The mean of the given variable subset.
    */
double calculate_mean(double** subset, const size_t subset_size);

/**
    * @brief Calculate the standard deviation of a given variable subset.
    * @param subset Variable observations.
    * @param X_mean Variable mean. 
    * @param subset_size Number of observations.
    * @return The standard deviation of the given variable subset.
    */
double calculate_standard_deviation(double** subset, const double subset_mean, const size_t subset_size);

double calculate_pearson_coeffient(double** X_subset, double** Y_subset, const size_t subset_syze)
{
    const double X_mean = calculate_mean(X_subset, subset_syze);

    const double Y_mean = calculate_mean(Y_subset, subset_syze);

    const double covariance = calculate_covariance(X_subset, X_mean, Y_subset, Y_mean, subset_syze);

    const double X_standard_deviation = calculate_standard_deviation(X_subset, X_mean, subset_syze);

    const double Y_standard_deviation = calculate_standard_deviation(Y_subset, Y_mean, subset_syze);

    return covariance / (X_standard_deviation * Y_standard_deviation) ;
}

double calculate_covariance(double** X_subset, const double X_mean, double** Y_subset, const double Y_mean, const size_t subset_syze)
{
    double sum = 0;

    for(size_t index = 0;  index < subset_syze; ++index)
    {
        sum += ( (*X_subset)[index] - X_mean ) * ( (*Y_subset)[index] - Y_mean );
    }

    return sum / (subset_syze - 1);
}

double calculate_mean(double** subset, const size_t subset_size)
{
    double sum = 0.0;

    for(size_t index = 0;  index < subset_size; ++index)
        sum += (*subset)[index];

    return sum / subset_size;
}

double calculate_standard_deviation(double** subset, const double subset_mean, const size_t subset_size)
{
    double sum = 0.0;

    for(size_t index = 0;  index < subset_size; ++index)
        sum += pow( (  (*subset)[index] - subset_mean  ), 2 );

    return sqrt( sum / (subset_size-1) );
}

bool is_correlated(const double pearson_correlation_coefficient, const double lower_bound, const double upper_bound)
{
    return pearson_correlation_coefficient >= lower_bound && pearson_correlation_coefficient <= upper_bound;
}
