#ifndef ARGS_H
#define ARGS_H

#include <stdbool.h>

typedef struct
{
	//~ //fields
	bool extended;
	bool ignore_case;
	bool new_line;
	const char* pattern;
	
	//mios
	bool transpose;
	bool corre;
	bool anti_corre;
	bool output;
	bool print;
	
	char *input_file;
	char *output_file;
	const char* cancer;
	char arguments[5][100];
	
}args_t;

/**
    * @brief Initializes the program arguments values.
    * @param args An struct where command line arguments will be stored.
    */

void args_init(args_t* args);

/**
    * @brief Validate the correctness of the command line arguments given by the user.
    * @param args An struct where command line arguments will be stored if no errors were detected.
    * @param argc Number of command line arguments.
    * @param argv Command line arguments.
    * @return EXIT_SUCESS if no errors were detected in the given command line arguments.
    */
int args_analyze(args_t* args, int argc, char ** argv);

/**
    * @brief Print the expected command line arguments format.
    * @return Zero to allow the user to restart the program execution with the correct parameters,
    */
int args_print_help(void);

/**
    * @brief Free the memory required to store the program arguments.
    * @param args Program arguments.
    */
void args_destroy(args_t* args);


#endif // ARGS_H
