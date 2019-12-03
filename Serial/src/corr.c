#include "corr.h"
#include "correlation_coefficient_summarizer.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/**
 * @brief Gets all the necessary flags that regular expressions need
 * @param corr Pointer to the class' struct 
 * @return cflags
 * */
 
int convert_cflags(const corr_t* corr);

/**
 * @brief Generates the output.csv file
 * @param corr Pointer to the class' struct 
 * @param correlation_record Array with the positions of all the correlated types of cancer
 * @param transpose True if the user used the -t flag
 * */
void generate_file(corr_t* corr, int* correlation_record, bool transpose);

/**
 * @brief Generates the output.csv file
 * @param corr Pointer to the class' struct 
 * @param matches Array used when user triggers the [regex] option.
 * */
void get_matches(corr_t* corr, int* matches);

/**
 * @brief Prints the matrix with all the Pearson's correlation coefficients
 * @param corr Pointer to the class' struct.
 * @param correlation_matrix Matrix to fill with all the correlation coefficients found when comparing
 * */
void print_correlation_matrix(corr_t* corr, double*** correlation_matrix);

void corr_init(corr_t* corr)
{
	args_init( &corr->args );
}

int corr_run(corr_t* corr, int argc, char ** argv)
{
	int* correlation_record;
	int* matches = NULL;
	
	int error = args_analyze( &corr->args, argc, argv );
	if( error )
		return error;
	
	if( corr->args.pattern == NULL)
		return args_print_help();	

	load_file(corr->args.input_file, &corr->data, corr->args.transpose);	
	correlation_record = (int*) calloc(corr->data.column_count, sizeof(int));
	
	double** correlation_coefficients = (double**) calloc( (corr->data.column_count-1) , sizeof(double*));
	for(int row = 0; row < (corr->data.column_count-1); ++row)
	{
		correlation_coefficients[row] = (double*) calloc( (corr->data.column_count-1), sizeof(double));
	}

	double upper_bound = 0.0;
	double lower_bound = 0.0;
		
	if( corr->args.corre )
	{
		if( strcmp(&corr->args.arguments[0][0],"0") ) 
		{
			upper_bound = atof(&corr->args.arguments[1][0]);
			lower_bound = atof(&corr->args.arguments[0][0]);
		}else
		{
			upper_bound = 1;
			lower_bound	= 0.75;
		}
			
	}else
	{	
		if(strcmp(&corr->args.arguments[2][0],"0"))	
		{
			upper_bound = atof(&corr->args.arguments[3][0]);
			lower_bound = atof(&corr->args.arguments[2][0]);
		}else
		{
			upper_bound = -1;
			lower_bound	= -0.75;
		}
	}
	
	if(corr->args.cancer){
		matches = (int*) calloc(corr->data.column_count, sizeof(int));
		int cflags = convert_cflags(corr);
	
		if( regcomp(&corr->regex , corr->args.cancer, cflags) )
			return fprintf(stderr, "error: invalid regular expression: %s\n", corr->args.cancer), 3;
		
		get_matches(corr, matches);		
	
		regfree( &corr->regex );
	}
		
	start_summarazing(&corr->data.values, &correlation_coefficients, corr->data.column_count-1, corr->data.row_count-1, lower_bound, upper_bound, &correlation_record, matches, corr->args.transpose);
	
	if(corr->args.print)
		print_correlation_matrix(corr, &correlation_coefficients);
	
	generate_file(corr, correlation_record, corr->args.transpose);
		
	
	for(int row = 0; row < (corr->data.column_count-1); ++row)
	{
		free(correlation_coefficients[row]);
	}
		
	free(correlation_coefficients);
	free(correlation_record);
	free(matches);
	corr_destroy(corr);
	return EXIT_SUCCESS;
}

void corr_destroy(corr_t* corr)
{
	args_destroy( &corr->args );
	parse_destroy( &corr->data, corr->args.transpose );
}

void get_matches(corr_t* corr, int* matches)
{
	int eflags = 0;
	
	for(int column = 0;column < corr->data.column_count-1; ++column)
	{
		if ( regexec( &corr->regex, &corr->data.names[column][0], 0, NULL, eflags) == 0 )
			matches[column] = 1;
	}
}


void print_correlation_matrix(corr_t* corr, double*** correlation_matrix)
{
	
	for(int row = 0; row < corr->data.column_count; ++row)
	{
		for(int column = 0;column < corr->data.column_count; ++column)
		{
			if( row == 0 && column ==  0)
				printf("[------],");
			else if( row == 0)
				printf("%s, ",&corr->data.names[column-1][0]);
			else if( column == 0 )
				printf("%s, ",&corr->data.names[row-1][0]);
			else
				printf("%lf, ", (*correlation_matrix)[row-1][column-1]);
				
		}
		printf("\n");
	}
}

void generate_file(corr_t* corr, int* correlation_record, bool transpose)
{
	FILE *file;
	char* write;
	
	file = fopen(corr->args.output_file, "w");
	
	for(int row = 0; row < corr->data.row_count; ++row)
	{
		for(int column = 0;column < corr->data.column_count; ++column)
		{
			if(row == 0 && column == 0)
			{
				fprintf(file, "%s,","");
			}else if(correlation_record[column] == 1)
			{
				if(row == 0)
				{
					if(!transpose)
						write = &corr->data.names[column-1][0];
					else
						write = &corr->data.gens[column-1][0];
				}else if(column == 0 && row != 0)
				{
					if(!transpose)
						write = &corr->data.gens[row-1][0];
					else
						write = &corr->data.names[row-1][0];
				}else
				{
					sprintf(write,"%f",corr->data.values[row-1][column-1]);
				}
				
				fprintf(file, "%s",write);
				
				if(column != corr->data.column_count-1)
					fprintf(file, ",");
			}

		}
		fprintf(file, "\n");
	}
	fclose(file);
}

int convert_cflags(const corr_t* corr)
{
	int cflags = REG_NOSUB;
	if( corr->args.extended ) 
		cflags |= REG_EXTENDED;
		
	if( corr->args.ignore_case ) 
		cflags |= REG_ICASE;

	if( corr->args.new_line ) 
		cflags |= REG_NEWLINE;
		
	return cflags;
}


