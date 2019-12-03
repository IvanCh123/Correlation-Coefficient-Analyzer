#ifndef CORR_H
#define CORR_H

#include "args.h"
#include "parse_file.h"
#include <regex.h>

typedef struct
{
	//fields
	args_t args;
	data_t data;
	regex_t regex;

}corr_t;

/**
    * @brief Initialize the program data, arguments, and regular expression match type.
    * @param corr An struct to store the program data, arguments, and regular expression match type.
    */
void corr_init(corr_t* corr);

/**
    * @brief It executes a program that receives a data set with the objective of reducing it using as 
    * a discard mechanism the correlation-anticorrelation between its variables.
    * @param corr An struct to store the program data, arguments, and regular expression match type.
    * @param argc Number of command line arguments.
    * @param argv Command line arguments.
    * @return EXIT_SUCESS if an error that prevented the correct execution of the program was not detected.
    */
int corr_run(corr_t* corr, int argc, char ** argv);

/**
    * @brief Free the memory required to store the program data, arguments and regular expression match type.
    * @param corr An struct containg the above program components.
    */
void corr_destroy(corr_t* corr);


#endif // CORR_H
