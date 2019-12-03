#include "corr.h"
 
int main(int argc, char ** argv){

	corr_t corr;
	corr_init(&corr);
	corr_run(&corr, argc, argv);
	
	return 0;	
}
