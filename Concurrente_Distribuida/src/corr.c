#include "corr.h"
#include "mathematical_operations.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <omp.h>
#include <mpi.h>

#define MIN(a,b) ( (a) < (b) ? (a) : (b))

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
 * */
void generate_file(corr_t* corr, int* correlation_record);

/**
 * @brief Generates the output.csv file
 * @param corr Pointer to the class' struct 
 * @param @param matches Array used when user triggers the [regex] option.
 * */
void get_matches(corr_t* corr, int* matches);

/**
 * @brief Prints the matrix with all the Pearson's correlation coefficients
 * @param corr Pointer to the class' struct.
 * @param correlation_matrix Matrix to fill with all the correlation coefficients found when comparing
 * */
void print_correlation_matrix(corr_t* corr, double*** correlation_matrix);

typedef struct 
{
    double lower_bound;         		// Correlation / anti-correlation lower bound.
    double upper_bound;         		// Correlation / anti-correlation upper bound.

}data_set_info_t;

/**
 * @brief Sets the upeer bound and the lower bound
 * @param corr Pointer to the class' struct.
 * @param info Pointer to the local struct.
 * */
void set_range(corr_t* corr, data_set_info_t* info);


/**
    * @brief Summarizes the given data set by checking the correlation coefficient between every pair of variables (cancer types) in the correlation matrix. 
    * If a variable is correlated-anticorrelated with at least one another it's conserved, otherwise discarded. We keep track of it my using an array
    * whose cells represent each variable and it'll be set to one if the corresponding columns must be conserved in the outputfile. 
    * @param corr Pointer to the class' struct.
    * @param correlation_coefficients Matrix to fill with all the correlation coefficients found when comparing
    * @param correlation_record Array to fill with 1's or 0's when a cancer type is correlated to another
    * @param start Column start
    * @param finish Column finish
    * @param my_rank Process ID
    */ 
void summarize_all_data(data_set_info_t* info, corr_t* corr, double*** correlation_coefficients, int* correlation_record, int start, int finish, int my_rank);


/**
    * @brief Summarizes the given data set by checking the correlation coefficient the specified cancer types in the correlation matrix. 
    * If a variable is correlated-anticorrelated with at least one another it's conserved, otherwise discarded. We keep track of it my using an array
    * whose cells represent each variable and it'll be set to one if the corresponding columns must be conserved in the outputfile. 
    * @param info A struct containing the data set to be summarized, it's dimensions.
    * @param corr Pointer to the class' struct.
    * @param correlation_coefficients Matrix to fill with all the correlation coefficients found when comparing
    * @param correlation_record Array to fill with 1's or 0's when a cancer type is correlated to another
    * @param matches Array used when user triggers the [regex] option.
    * @param start Column start
    * @param finish Column finish
    */ 
void summarize_specified_data(data_set_info_t* info, corr_t* corr, double*** correlation_coefficients, int* correlation_record, int* matches, int start, int finish, int my_rank);

/**
    * @brief Calculates and stores the correlation coefficient between every pair of variables (cancer types). 
    * @param corr Pointer to the class' struct.
    * @param values Matrix with all the floating points values the input file has
    * @param correlation_coefficients Matrix to fill with all the correlation coefficients found when comparing
    * @param start Column start
    * @param finish Column finish
    * @param my_rank Process ID
    * @param process_count Amounts of processors
    */
void fill_correlation_matrix(corr_t *corr, double** values, double*** correlation_coefficients, const int start, const int finish, int my_rank, int process_count);

/***
 * @brief Sends all the values stored in the matrix of correlation coefficients
 * @param corr Pointer to the class' struct.
 * @param correlation_coefficients Matrix to fill with all the correlation coefficients found when comparing	
 * @param process_count Amounts of processors
 * @param my_rank Process ID
 * */
void share_matrix(corr_t* corr, double*** correlation_coefficients, int process_count, int my_rank);

/**
    * @brief Summarize the given data subset using as a discard mechanism the correlation or anticorrelation between its variables (cancer types).
    * @param corr Pointer to the class' struct.
    * @param values Matrix with all the floating points values the input file has
    * @param correlation_coefficients Matrix to fill with all the correlation coefficients found when comparing
    * @param correlation_record Array to fill with 1's or 0's when a cancer type is correlated to another
    * @param matches Array used when user triggers the [regex] option. 
    * @param start Column start
    * @param finish Column finish
    * @param my_rank Process ID
    * @param process_count Amounts of processors
    */
void start_summarazing(data_set_info_t* info, corr_t *corr, double** values, double*** correlation_coefficients, int* correlation_record, int* matches, const int start, const int finish, int my_rank, int process_count );


int calculate_start( int data_count, int process_count, int process_id )
{
	return process_id * (data_count / process_count) + MIN( process_id, data_count % process_count);
}

int calculate_finish( int data_count, int process_count, int process_id )
{
	return calculate_start(data_count, process_count, process_id+1);
}


void corr_init(corr_t* corr)
{
	args_init( &corr->args );
}

int corr_run(corr_t* corr, int argc, char ** argv)
{
	MPI_Init(&argc, &argv);
	
	data_set_info_t info;
        
	int my_rank = -1;
	int process_count = -1;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_count);
	
	int error = args_analyze( &corr->args, argc, argv );
	
	if( error )
		return error;
		
	if( corr->args.pattern == NULL){
		if(my_rank == 0)
			return args_print_help();	
		else
			return 0;
	}
	
	load_file(corr->args.input_file, &corr->csv, corr->args.transpose);	
	
	if(my_rank==0)
		printf("Reading file: %s \n",corr->args.input_file);
	
	set_range(corr,&info);
			
	int* correlation_record = (int*) calloc(corr->csv.column_count, sizeof(int));
	
	double** correlation_coefficients = (double**) calloc( (corr->csv.column_count-1) , sizeof(double*));
	for(int row = 0; row < (corr->csv.column_count-1); ++row)
	{
		correlation_coefficients[row] = (double*) calloc( (corr->csv.column_count-1), sizeof(double));
	}
	
	int* matches = NULL;
	if(corr->args.cancer != NULL){
		matches = (int*) calloc(corr->csv.column_count, sizeof(int));
		int cflags = convert_cflags(corr);
		if( regcomp(&corr->regex , corr->args.cancer, cflags) )
		return fprintf(stderr, "error: invalid regular expression: %s\n", corr->args.cancer), 3;
		
		get_matches(corr, matches);		
	
		regfree( &corr->regex );
	}	
	
	int start = calculate_start(corr->csv.column_count-1, process_count, my_rank);
	
	int finish = calculate_finish(corr->csv.column_count-1, process_count, my_rank);
	
	
	if(corr->csv.column_count < process_count)
		process_count = corr->csv.column_count;
	

	
	
	start_summarazing(&info, corr, corr->csv.values, &correlation_coefficients, correlation_record, matches, start, finish, my_rank, process_count);
	

	if(my_rank == 0 && corr->args.print)
		print_correlation_matrix(corr, &correlation_coefficients);
	
	if(my_rank == 0)
		generate_file(corr, correlation_record);
	
	

	free(matches);
	for(int pos=0;pos<corr->csv.column_count-1;++pos)
		free(correlation_coefficients[pos]);	
	free(correlation_coefficients);	
	free(correlation_record);
	corr_destroy(corr);
	
	MPI_Finalize();	
	
	return EXIT_SUCCESS;
}

void start_summarazing(data_set_info_t* info, corr_t* corr, double** values, double*** correlation_coefficients, int* correlation_record, int* matches, const int start, const int finish, int my_rank, int process_count )
{
	fill_correlation_matrix(corr, values, correlation_coefficients, start, finish, my_rank, process_count);
		
	share_matrix(corr, correlation_coefficients, process_count, my_rank);
	
	if( matches == NULL )
		summarize_all_data(info, corr, correlation_coefficients, correlation_record, start,finish,my_rank);
	else
		summarize_specified_data(info, corr, correlation_coefficients, correlation_record, matches, start, finish, my_rank);
}

void fill_correlation_matrix(corr_t *corr, double** values, double*** correlation_coefficients, const int start, const int finish, int my_rank, int process_count )
{
	// The pearson's correlation coefficient between every variables X and every variable Y is calculated.
    // If a variable is not correlated with any other it is discarded. For this, a record is kept in an array that
    // indicates with a 1 if the variable corresponding to its index was correlated with at least one variable and a 0 if not.

	double* X_variable_subset = (double*) calloc( corr->csv.row_count-1, sizeof(double) );  // Temporal array for X_subset.
	double* Y_variable_subset = (double*) calloc( corr->csv.row_count-1, sizeof(double) );  // Temporal array for Y_subset.
	
	#pragma omp parallel for num_threads(process_count) default(none) shared(X_variable_subset, Y_variable_subset,process_count, corr, values, correlation_coefficients, my_rank)
	for(int X_variable = start; X_variable < finish; ++X_variable)  // First cancer type.
	{		
		for(int Y_variable = 0; Y_variable < corr->csv.column_count-1; ++Y_variable) // Second cancer type.
		{	
			for(int observation = 0; observation < corr->csv.row_count-1; ++observation) // Current gen.
			{
				X_variable_subset[observation] =  values[observation][X_variable];
				Y_variable_subset[observation] =  values[observation][Y_variable];
			}
			
			double val = 0.0;
			val = calculate_pearson_coeffient( &X_variable_subset, &Y_variable_subset, corr->csv.row_count-1 );
			
			if(my_rank != 0){
				MPI_Send(&X_variable,1,MPI_INT,0,0,MPI_COMM_WORLD);
				MPI_Send(&Y_variable,1,MPI_INT,0,0,MPI_COMM_WORLD);
				MPI_Send(&val,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
			}else{
				int posX =0 ;
				int posY = 0;
				double val_r=0.0;
					
				(*correlation_coefficients)[X_variable][Y_variable] = val;	
				for(int i=1;i<process_count;i++){
					MPI_Recv(&posX,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					MPI_Recv(&posY,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					MPI_Recv(&val_r,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					(*correlation_coefficients)[posX][posY] = val_r;
			    }	    
			}	
				
		}
				
	}
	
	free(X_variable_subset);
	free(Y_variable_subset);
}

void share_matrix(corr_t* corr, double*** correlation_coefficients, int process_count, int my_rank)
{
	MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank == 0){
		double val_s = 0.0;
		
		for(int row=0;row<corr->csv.column_count-1;++row)
		{
			for(int column=0;column<corr->csv.column_count-1;++column)
			{
				for(int process = 1;process<process_count;++process)
				{	
					val_s = (*correlation_coefficients)[row][column];
					MPI_Send(&row,1,MPI_INT,process,0,MPI_COMM_WORLD);
					MPI_Send(&column,1,MPI_INT,process,0,MPI_COMM_WORLD);
					MPI_Send(&val_s,1,MPI_DOUBLE,process,0,MPI_COMM_WORLD);
				}
			}
		}
		
	}else{
		double val_r = 0.0;
		int posX = 0;
		int posY = 0;
		
			for(int row=0;row<corr->csv.column_count-1;++row)
			{
				for(int column=0;column<corr->csv.column_count-1;++column)
				{

					MPI_Recv(&posX,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					MPI_Recv(&posY,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					MPI_Recv(&val_r,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					(*correlation_coefficients)[posX][posY] = val_r;
				}
			}
	}
}

void summarize_all_data(data_set_info_t* info, corr_t* corr, double*** correlation_coefficients, int* correlation_record, int start, int finish, int my_rank)
{
	correlation_record[0] = 1;
	
	int X_variable;
	int Y_variable;
	
	#pragma omp parallel for default(none) shared(info,correlation_record,corr,correlation_coefficients, start, finish, my_rank, X_variable, Y_variable)
	for(X_variable = start; X_variable < finish; ++X_variable)  // First cancer type.
	{
		for(Y_variable = X_variable+1; Y_variable < corr->csv.column_count-1; ++Y_variable) // Second cancer type.
		{
			double val = (*correlation_coefficients)[X_variable][Y_variable];
			if( is_correlated(val, info->lower_bound, info->upper_bound) )
			{
					
				if(my_rank != 0){
					int pos1= X_variable+1;
					int pos2= Y_variable+1;
					MPI_Send(&pos1,1,MPI_INT,0,0,MPI_COMM_WORLD);
					MPI_Send(&pos2,1,MPI_INT,0,0,MPI_COMM_WORLD);
					
				}else{
					
					int posX =X_variable ;
					int posY = Y_variable;
					#pragma omp critical
					{
						correlation_record[posX] = 1;
						correlation_record[posY] = 1;
					}
					#pragma omp critical
					{
						int p1 = 0;
						int p2 = 0;
						MPI_Recv(&p1,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
						MPI_Recv(&p2,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

						correlation_record[p1] = 1;	
						correlation_record[p2] = 1;
					}	
				}	
			}
		}
	}
}


void summarize_specified_data(data_set_info_t* info, corr_t* corr, double*** correlation_coefficients, int* correlation_record, int* matches, int start, int finish, int my_rank)
{
	correlation_record[0] = 1;
	
	int X_variable;
	int Y_variable;
	#pragma omp parallel for default(none) shared(info,correlation_record,corr,correlation_coefficients, matches, start, finish, my_rank, X_variable, Y_variable) 
	for(X_variable = start; X_variable < finish; ++X_variable)  // First cancer type.
	{
		if( matches[X_variable] )
		{
			for(Y_variable = 0; Y_variable < corr->csv.column_count-1; ++Y_variable) // Second cancer type.
			{
				if( X_variable != Y_variable )
				{
					double val = (*correlation_coefficients)[X_variable][Y_variable];
					if( is_correlated(val, info->lower_bound, info->upper_bound ) )
					{
						if(my_rank != 0){
							int pos1= X_variable+1;
							int pos2= Y_variable+1;
							MPI_Send(&pos1,1,MPI_INT,0,0,MPI_COMM_WORLD);
							MPI_Send(&pos2,1,MPI_INT,0,0,MPI_COMM_WORLD);
						}else{
					
							int posX =X_variable ;
							int posY = Y_variable;
							#pragma omp critical
							{
								correlation_record[posX] = 1;
								correlation_record[posY] = 1;
							}
							#pragma omp critical
							{
								int p1 = 0;
								int p2 = 0;
								MPI_Recv(&p1,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
								MPI_Recv(&p2,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

								correlation_record[p1] = 1;	
								correlation_record[p2] = 1;
							}	
						}	
					}
				}
			}
		}
	}
}

void corr_destroy(corr_t* corr)
{
	args_destroy( &corr->args );
	parse_destroy( &corr->csv, corr->args.transpose );
}

void get_matches(corr_t* corr, int* matches)
{
	int eflags = 0;
	
	for(int column = 0;column < corr->csv.column_count-1; ++column)
	{
		if ( regexec( &corr->regex, &corr->csv.names[column][0], 0, NULL, eflags) == 0 )
			matches[column] = 1;
	}
}

void set_range(corr_t* corr,data_set_info_t* info)
{
	if( corr->args.corre )
	{	
		if( strcmp(&corr->args.arguments[0][0],"0") ) 
		{
			info->upper_bound = atof(&corr->args.arguments[1][0]);
			info->lower_bound = atof(&corr->args.arguments[0][0]);
		}else
		{
			info->upper_bound = 1;
			info->lower_bound	= 0.75;
		}	
	}else
	{	
		if(strcmp(&corr->args.arguments[2][0],"0"))	
		{
			info->upper_bound = atof(&corr->args.arguments[3][0]);
			info->lower_bound = atof(&corr->args.arguments[2][0]);
		}else
		{
			info->upper_bound = -1;
			info->lower_bound = -0.75;
		}
	}	
}


void print_correlation_matrix(corr_t* corr, double*** correlation_matrix)
{
	for(int row = 0; row < corr->csv.column_count; ++row)
	{
		for(int column = 0;column < corr->csv.column_count; ++column)
		{
			if( row == 0 && column ==  0)
				printf("[------],");
			else if( row == 0)
				printf("%s, ",&corr->csv.names[column-1][0]);
			else if( column == 0 )
				printf("%s, ",&corr->csv.names[row-1][0]);
			else
				printf("%lf, ", (*correlation_matrix)[row-1][column-1]);
				
		}
		printf("\n");
	}
}

void generate_file(corr_t* corr, int* correlation_record)
{
	FILE *file;
	char* write;
	
	file = fopen(corr->args.output_file, "w");
	
	for(int row = 0; row < corr->csv.row_count; ++row)
	{
		for(int column = 0;column < corr->csv.column_count; ++column)
		{
			if(row == 0 && column == 0)
			{
				fprintf(file, "%s,","");
			}else if(correlation_record[column] == 1)
			{
				if(row == 0)
				{
					write = &corr->csv.names[column-1][0];

				}else if(column == 0 && row != 0)
				{
					write = &corr->csv.gens[row-1][0];
				}else
				{
					sprintf(write,"%f",corr->csv.values[row-1][column-1]);
				}
				
				fprintf(file, "%s",write);
				
				if(column != corr->csv.column_count-1)
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


