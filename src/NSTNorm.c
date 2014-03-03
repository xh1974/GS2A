/*
 *  NSTNorm.c
 *	Normal Score Transform applied to expression data
 *
 *  Created by Han Xu on 02/12/13.
 *  Copyright 2013 Dana Farber Cancer Institute. All rights reserved.
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include "math_api.h"
#include "dataMatrix.h"


//Normalize data using Normal Score Transformation
void NSTNormData(double *data, int geneNum, int sampleNum);

//print command usage 
void PrintCommandUsage();


//Normalize data using Normal Score Transformation
void NSTNormData(double *data, int geneNum, int sampleNum)
{
	int i;
	int *rank;
	double *normData;
	
	rank = (int *)malloc(sampleNum*sizeof(double));
	normData = (double *)malloc(sampleNum*sizeof(double));
	
	assert((rank!=NULL)&&(normData!=NULL));
							
	for (i=0;i<geneNum;i++)
	{
		Ranking(rank, data+i*sampleNum, sampleNum);
		NormalTransform(normData, rank, sampleNum);
		memcpy(data+i*sampleNum, normData, sampleNum*sizeof(double));
	}
	
	free(rank);
	free(normData);
}

//print command usage 
void PrintCommandUsage()
{
	//print the options of the command
	printf("NSTNorm - Normal Score Transform for expression data.\n");
	printf("usage:\n");
	printf("-i <input data matrix>\n");
	printf("-o <output normalized data matrix>\n");
	printf("example:\n");
	printf("NSTNorm -i input.txt -o output.txt \n");
}

int main (int argc, const char * argv[]) 
{
	char expressionFileName[1000], outputFileName[1000];
	DATA_MATRIX_STRUCT expressions;
	int i;
	
	
	//Parse the command line
	if (argc == 1)
	{
		PrintCommandUsage();
		return -1;
	}
	
	expressionFileName[0] = 0;
	outputFileName[0] = 0;
	
	for (i=2;i<argc;i++)
	{
		if (strcmp(argv[i-1], "-i")==0)
		{
			strcpy(expressionFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-o")==0)
		{
			strcpy(outputFileName, argv[i]);
		}
	}
	
	if ((expressionFileName[0]==0)||(outputFileName[0]==0))
	{
		printf("Command error!\n");
		PrintCommandUsage();
		return -1;
	}
	
	if (ReadDataMatrix(expressionFileName, &expressions)<=0)
	{
		printf("Cannot open %s or file format error!\n", expressionFileName);
		return -1;
	}
	
	printf("sampleNum=%d\ngeneNum=%d\n", expressions.sampleNum, expressions.recordNum);
	
	NSTNormData(expressions.matrix, expressions.recordNum, expressions.sampleNum);
	
	if (SaveDataMatrix(outputFileName, &expressions)<=0)
	{
		printf("Cannot save to %s\n",outputFileName);
	}
	else
	{
		printf("NST finished.\n");
	}
	
	FreeDataMatrix(&expressions);
	return 0;
}
