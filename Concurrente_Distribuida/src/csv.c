#include "csv.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NAME_CAPACITY 1000
#define LINE_CAPACITY 10000

void load_file(char *input_file, csv_t* data, bool transpose)
{	
	FILE* file;
	
	char singleLine[LINE_CAPACITY];

	int row=0;
	int column=0;
	int character_count=0;

	file = fopen(input_file, "r");
	
	while(!feof(file))
	{		
		character_count= fgetc(file);
		if(character_count == ',')
			++column;
			
		if(character_count=='\n')	
		{
			data->column_count = column+1;
			column = 0;
			++row;
		}
			
	}
	fclose(file);
	
	data->row_count=row;
	data->names = (char**) calloc (data->column_count-1, sizeof(char*));
	data->gens = (char**) calloc (data->row_count, sizeof(char*));
	
	if( transpose )
	{
		int temp = data->row_count;
		data->row_count = data->column_count;
		data->column_count = temp;
		
		data->values = (double**) calloc (data->row_count, sizeof(double*));
		
		for(int index = 0; index < data->row_count; ++index)
		{
			data->values[index] = (double*) calloc (data->column_count, sizeof(double));	
			data->names[index] = (char*) calloc (NAME_CAPACITY, sizeof(char));	
			
		}
		for(int index = 0; index < data->column_count; ++index)
		{
			data->gens[index] = (char*) calloc (NAME_CAPACITY, sizeof(char));
		}
		
	}else{
		data->values = (double**) calloc (data->row_count, sizeof(double*));
		
		for(int index = 0; index < data->row_count; ++index)
		{
			data->values[index] = (double*) calloc (data->column_count, sizeof(double));		
			data->gens[index] = (char*) calloc (NAME_CAPACITY, sizeof(char));
		}
		
		for(int index = 0; index < data->column_count-1; ++index)
		{
			data->names[index] = (char*) calloc (NAME_CAPACITY, sizeof(char));
		}
	}
		
	row=0;
	column=0;
	int posNames=0;
	int posGens=0;
	int posArray=0;
	
	char* token;
	file = fopen(input_file, "r");
	while(!feof(file) && row <= data->row_count)
	{
		fgets(singleLine, LINE_CAPACITY, file);
		token = strtok(singleLine, ",");
		
		while( token != NULL) 
		{	
			if(row == 0 )
			{
				int si = strlen(token) - 1;
				if(token[si] == '\n')
					token[si] = '\0';
				strcpy(&data->names[posNames][0], token);
				++posNames;	
					
			}
			else
			{
				if(posArray == 0)
				{
					strcpy(&data->gens[posGens++][0], token);				
				}
				else
				{
					if( !transpose )
						data->values[row-1][column++] = atof(token);	
					else
						data->values[column++][row-1] = atof(token);	
				}

			}
			token = strtok(NULL, ",");
			++posArray;
		}
		posArray=0;
		column=0;
		++row;
		
	}
	fclose(file);
	
}

void parse_destroy(csv_t* data, bool transpose){
	
	if(!transpose)
	{
		for(int index = 0; index<data->row_count;++index){
			free( data->values[index] );
			free( data->gens[index] );
		}	
		for(int index = 0; index<data->column_count;++index)
			free( data->names[index] );
		
	}else{
		
		for(int index = 0; index<data->row_count;++index){
			free( data->names[index] );
			free( data->values[index] );
		}
		for(int index = 0; index<data->column_count;++index){
			free( data->gens[index] );
			
		}
	}

	free(data->gens);
	free(data->names);
	free(data->values);
}
