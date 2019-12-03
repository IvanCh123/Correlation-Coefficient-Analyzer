#include "args.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char* const help = 
	"expected command: ./corr input_file.csv -cc/ac [range] [opts] -o output_file.csv [regex]\n"
	"options : \n"
	"	-e  extended\n"
	"	-i  ignore case\n"
	"	-n  new line chars separate strings\n"
	"   -m  Print correlation matrix to FILE\n"
	"	-t  transpose output matrix\n"
	"	-cc correlation\n"
	"	-ac anti correlation\n"
	"	-o  [output_file] output file\n"
	;

void args_init(args_t* args)
{	
	args->extended = false;
	args->ignore_case = false;
	args->new_line = false;
	
	
	args->transpose = false;
	args->corre = false;
	args->anti_corre = false;
	args->output = false;
	args->print = false;
	
	args->pattern = NULL;
	args->input_file = NULL;
	args->output_file = NULL;
	args->cancer = NULL;
}

int args_analyze(args_t* args, int argc, char ** argv)
{
	char convert[100];
	char *token1;
	double upper_r = 0.0;
	double lower_r = 0.0;
	
	for(int index = 1; index < argc; ++index)
	{
		if( *argv[index] == '-' ) //si el argumento empieza con - es opcion
		{
			for ( const char* option = argv[index] + 1; *option; ++option)
			{
				switch( argv[index][1] )
				{	
					case 'e': args->extended = true; break;
					
					case 'i': args->ignore_case = true; break;
					
					case 'n': args->new_line = true; break;	
					
					case 'm': args->print = true; break;	
							
					case 't': args->transpose = true; break;	
									
					case 'c': 
						if(!args->anti_corre)
						{
							if(argv[index+1] != NULL && *argv[index+1] != '-')
							{
								sscanf(argv[index+1],"%lf:%lf",&upper_r,&lower_r);
								if(upper_r >= 0.0 && lower_r >= 0.0)
								{
									sprintf(convert,"%lf",upper_r);
									strcpy(&args->arguments[0][0], convert);
									sprintf(convert,"%lf",lower_r);
									strcpy(&args->arguments[1][0], convert);
								}else{
									return fprintf(stderr, "error: The correlation range must be positive\n"), EXIT_FAILURE; 	
								}	
								
								
							}else{
								if(argv[index+1] != NULL && strcmp(argv[index+1],"-o"))
									return fprintf(stderr, "error: The correlation range must be positive\n"), EXIT_FAILURE; 	
								strcpy(&args->arguments[0][0], "0");
							}
							
							args->corre = true; 
						}else{
							return fprintf(stderr, "error: You can't choose -cc and -ac at the same time."), EXIT_FAILURE; 
						}
						break;	
									
					case 'a': 
						if(!args->corre)
						{
							if( argv[index+1] != NULL && *argv[index+1] == '-' && strcmp(argv[index+1],"-o") )
							{
								sscanf(argv[index+1],"%lf:%lf",&upper_r,&lower_r);

								if(upper_r <= 0.0 && lower_r <= 0.0)
								{
									sprintf(convert,"%lf",upper_r);
									strcpy(&args->arguments[2][0], convert);
									sprintf(convert,"%lf",lower_r);
									strcpy(&args->arguments[3][0], convert);
									index+=2;
								}else{
									return fprintf(stderr, "error:1 The anti-correlation range must be negative\n"), EXIT_FAILURE; 	
								}	
							}else{
								if(argv[index+1] != NULL && *argv[index+1] != '-')
									return fprintf(stderr, "error:1 The anti-correlation range must be negative\n"), EXIT_FAILURE; 	
								strcpy(&args->arguments[2][0], "0");
							}
							
							
							args->anti_corre = true; 
						}else
						{
							return fprintf(stderr, "error: You can't choose -cc and -ac at the same time."), EXIT_FAILURE; 
						}
						break;	
									
					case 'o': 
						if(argv[index+1]==NULL)
							return fprintf(stderr, "You must indicate the output file with the following syntaxis: -o [output_file.csv]\n"), EXIT_FAILURE;	
						
						strcpy(convert,argv[index+1]);
						token1 = strtok(convert, ".");
						token1 = strtok(NULL, ".");
						if(token1 != NULL){
							if( !strcmp(token1,"csv") )
							{
								args->output_file = argv[index+1];
								args->output = true; 
							}else
								return fprintf(stderr, "You must indicate the output file with the following syntaxis: -o [output_file.csv]\n"), EXIT_FAILURE;	
						}
						break;				
					
					default: return fprintf(stderr, "error: invalid option: %c\n", *option), EXIT_FAILURE;
				}
			}
		}
		else
		{
			if( (args->corre || args->anti_corre) && args->output )
			{
				args->pattern = argv[0];
				if(index == argc-1 && strcmp(args->output_file,argv[index]))
					args->cancer = argv[index];
			}
			
			if(index == 1)
				args->input_file =  argv[index];
		}
	}
	if( !args->output && args->pattern )
		return fprintf(stderr, "You must indicate the output file with the following syntaxis: -o [output_file.csv]\n"), EXIT_FAILURE;			
	
	return EXIT_SUCCESS;
}

void args_destroy(args_t* args)
{
	(void)args;
}


int args_print_help(void)
{
	printf("%s", help);
	return 0;
}
