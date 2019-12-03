#include "correlation_coefficient_summarizer.h"

#include <stdlib.h>

/**
    * @brief Summarizes the given data set by checking the correlation coefficient between every pair of variables (cancer types) in the correlation matrix. 
    * If a variable is correlated-anticorrelated with at least one another it's conserved, otherwise discarded. We keep track of it my using an array
    * whose cells represent each variable and it'll be set to one if the corresponding columns must be conserved in the outputfile. 
    * @param info A struct containing the data set to be summarized, it's dimensions.
    */ 
void summarize_all_data(data_set_info_t* info);

/**
    * @brief Summarizes the given data set by checking the correlation coefficient the specified cancer types in the correlation matrix. 
    * If a variable is correlated-anticorrelated with at least one another it's conserved, otherwise discarded. We keep track of it my using an array
    * whose cells represent each variable and it'll be set to one if the corresponding columns must be conserved in the outputfile. 
    * @param info A struct containing the data set to be summarized, it's dimensions.
    */ 
void summarize_specified_data(data_set_info_t* info);

/**
    * @brief Calculates and stores the correlation coefficient between every pair of variables (cancer types). 
    * @param info A struct containing the correlation matrix to fill and it's dimensions.
    */
void fill_correlation_matrix(data_set_info_t* info);


/**
    * @brief Initialize the summarizer with the values given by the user.
    * @param data_set Complete data set to be reduced.
    * @param correlation_coeficients A matrix whose rows and columns associate two variables (cancer types) with their corresponding correlation coefficient.
    * @param variable_count Number of variables (cancer types).
    * @param subset_size Number of observations for each variable (gen types).
    * @param lower_bound Lower bound to determinate whether two variables are correlated-anticorrelated or not. 
    * @param upper_bound Upper bound to determinate whether two variables are correlated-anticorrelated or not. 
    * @param correlation_record An array whose cells correspond to each cancer type. It'll be one if the corresponding cancer type must be conserved, zero otherwise.
    * @param matches An array whose cells correspond to each cancer type. It'll be one if the corresponding cancer type matches the wildcard given by the user.
    * @param data_set_transposed True if the given data set is transposed.
    */
    
data_set_info_t get_data_set_info(double*** data_set, double*** correlation_coefficients, const size_t variable_count, const size_t subset_size, const double lower_bound, const double upper_bound, int** correlation_record, int* matches, const bool data_set_transposed);

struct data_set_info_t
{
    double** data_set;          		// Complete data set.
	double** correlation_coefficients;	// Correlation matrix.
    size_t variable_count;      		// Cancer type count (columns).
    size_t subset_size;         		// Gen type count (rows).
    double lower_bound;         		// Correlation / anti-correlation lower bound.
    double upper_bound;         		// Correlation / anti-correlation upper bound.
    int* correlation_record;    		// Each cell corresponds to a type of cancer. It is 1 if this type of cancer is correlated with at least one other and 0 if not.
	int* matches;						// For specified regular expressions.
	bool data_set_transposed;			// Transposed data set or not.
};

void start_summarazing(double*** data_set, double*** correlation_coefficients,const size_t variable_count, const size_t subset_size, const double lower_bound, const double upper_bound, int** correlation_record, int* matches,const bool data_set_transposed)
{
    data_set_info_t  info = get_data_set_info(data_set, correlation_coefficients, variable_count, subset_size, lower_bound, upper_bound,correlation_record, matches, data_set_transposed);	
	fill_correlation_matrix(&info);
	
	if( info.matches == NULL )
		summarize_all_data(&info);
	else
		summarize_specified_data(&info);
}

data_set_info_t get_data_set_info(double*** data_set, double*** correlation_coefficients , const size_t variable_count, const size_t subset_size, const double lower_bound, const double upper_bound, int** correlation_record, int* matches, const bool data_set_transposed)
{
    data_set_info_t info;
    info.data_set = *data_set;
    info.correlation_coefficients =  *correlation_coefficients;
    info.variable_count = (data_set_transposed) ? subset_size:variable_count;
    info.subset_size = (data_set_transposed) ? variable_count:subset_size;
    info.lower_bound = lower_bound;
    info.upper_bound = upper_bound;
    info.correlation_record = *correlation_record;
    info.matches = matches;
    info.data_set_transposed = data_set_transposed;
    return info;
}


void summarize_all_data(data_set_info_t* info)
{
	info->correlation_record[0] = 1;
	
	for(size_t X_variable = 0; X_variable < info->variable_count; ++X_variable)  // First cancer type.
	{
		for(size_t Y_variable = X_variable+1; Y_variable < info->variable_count; ++Y_variable) // Second cancer type.
		{
			if( is_correlated(info->correlation_coefficients[X_variable][Y_variable], info->lower_bound, info->upper_bound) )
			{
				info->correlation_record[X_variable+1] = 1;
				info->correlation_record[Y_variable+1] = 1;
			}
		}
	}
}

void summarize_specified_data(data_set_info_t* info)
{
	info->correlation_record[0] = 1;
	
	for(size_t X_variable = 0; X_variable < info->variable_count; ++X_variable)  // First cancer type.
	{
		if( info->matches[X_variable] )
		{
			for(size_t Y_variable = 0; Y_variable < info->variable_count; ++Y_variable) // Second cancer type.
			{
				if( X_variable != Y_variable )
				{
					if( is_correlated(info->correlation_coefficients[X_variable][Y_variable], info->lower_bound, info->upper_bound ) )
					{
						info->correlation_record[X_variable+1] = 1;
						info->correlation_record[Y_variable+1] = 1;
					}
				}
			}
		}
	}
}


void fill_correlation_matrix(data_set_info_t* info)
{
	double correlation_coefficient = 0.0;
	
	// The pearson's correlation coefficient between every variables X and every variable Y is calculated.
    // If a variable is not correlated with any other it is discarded. For this, a record is kept in an array that
    // indicates with a 1 if the variable corresponding to its index was correlated with at least one variable and a 0 if not.

	if( info->data_set_transposed )
	{
		for(size_t X_variable = 0; X_variable < info->variable_count; ++X_variable)  // First cancer type.
		{
			for(size_t Y_variable = 0; Y_variable < info->variable_count; ++Y_variable) // Second cancer type.
			{
				correlation_coefficient = calculate_pearson_coeffient( &info->data_set[X_variable], &info->data_set[Y_variable], info->subset_size );
				info->correlation_coefficients[X_variable][Y_variable] = correlation_coefficient;
				info->correlation_coefficients[Y_variable][X_variable] = correlation_coefficient;
			}
		}
	}	
	else
	{
		double* X_variable_subset = (double*) calloc( info->subset_size, sizeof(double) );  // Temporal array for X_subset.
		double* Y_variable_subset = (double*) calloc( info->subset_size, sizeof(double) );  // Temporal array for Y_subset.
		
		for(size_t X_variable = 0; X_variable < info->variable_count; ++X_variable)  // First cancer type.
		{
			for(size_t Y_variable = 0; Y_variable < info->variable_count; ++Y_variable) // Second cancer type.
			{
				for(size_t observation = 0; observation < info->subset_size; ++observation) // Current gen.
				{
					X_variable_subset[observation] =  info->data_set[observation][X_variable];
					Y_variable_subset[observation] =  info->data_set[observation][Y_variable];
				}
				
				correlation_coefficient = calculate_pearson_coeffient( &X_variable_subset, &Y_variable_subset, info->subset_size );
				info->correlation_coefficients[X_variable][Y_variable] = correlation_coefficient;
				info->correlation_coefficients[Y_variable][X_variable] = correlation_coefficient;
			
			}
		}
		free(X_variable_subset);
		free(Y_variable_subset);
	}
}
