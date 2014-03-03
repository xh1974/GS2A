/*
 *  dataMatrix.c
 *  processing data matrix
 *
 *  Created by Han Xu on 4/9/13.
 *  Copyright 2013 Dana Farber Cancer Institute. All rights reserved.
 *
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "dataMatrix.h"
#include "words.h"

//allocate memory for data matrix
int AllocDataMatrix(DATA_MATRIX_STRUCT *matrix, int sampleNum, int recordNum)
{
	matrix->recordInfo = (ID_INFO_STRUCT *)malloc(recordNum*sizeof(ID_INFO_STRUCT));
	matrix->sampleInfo = (ID_INFO_STRUCT *)malloc(sampleNum*sizeof(ID_INFO_STRUCT));	
	matrix->matrix = (double *)malloc(sampleNum*recordNum*sizeof(double));
	matrix->sampleNum = sampleNum;
	matrix->recordNum = recordNum;
	
	assert((matrix->sampleInfo!=NULL)&&(matrix->recordInfo!=NULL)&&(matrix->matrix!=NULL));
	
	if ((matrix->sampleInfo==NULL)||(matrix->recordInfo==NULL)||(matrix->matrix==NULL))
	{
		return -1;
	}
	else
	{
		return 1;
	}
}

//Free memory for data matrix
void FreeDataMatrix(DATA_MATRIX_STRUCT *matrix)
{
	free(matrix->sampleInfo);
	free(matrix->recordInfo);
	free(matrix->matrix);
	
	matrix->sampleNum = 0;
	matrix->recordNum = 0;
}


//Read the data matrix from a file
int ReadDataMatrix(char *fileName, DATA_MATRIX_STRUCT *matrix)
{
	FILE *fh;
	char **words;
	int wordNum;
	char *tmpS;
	int sampleNum, recordNum;
	int i;
	
	fh = (FILE *)fopen(fileName, "r");
	
	assert(fh!=NULL);
	
	if (!fh)
	{
		return -1;
	}
	
	words = AllocWords(MAX_SAMPLE_NUM+1, MAX_WORD_SIZE+1);
	
	assert(words!=NULL);
	
	if (words == NULL)
	{
		return -1;
	}
	
	tmpS = (char *)malloc((MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char));
	
	assert(tmpS!=NULL);
	
	recordNum = 0;
	
	//Read the header row to get the sample number
	fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
	
	wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");
	
	sampleNum = wordNum-1;
	
	//Read through the file to get the number of records
	
	fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
	wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");
	
	while ((wordNum==sampleNum+1)&&(!feof(fh)))
	{
		recordNum++;
		
		fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
		wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");
	}
	
	fclose(fh);
	
//	if ((recordNum<=0)||(sampleNum<=0)||(AllocDataMatrix(matrix, sampleNum, recordNum)<=0))
	if (AllocDataMatrix(matrix, sampleNum, recordNum)<=0)
	{		
		free(tmpS);
		FreeWords(words, MAX_SAMPLE_NUM+1);
		return -1;
	}
	
	//read the file again to retrieve the values
	
	fh = (FILE *)fopen(fileName, "r");
	
	assert(fh!=NULL);
	
	fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
	wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");

	for (i=0;i<sampleNum;i++)
	{
		strcpy(matrix->sampleInfo[i].name, words[i+1]);
	}
	
	recordNum = 0;
	
	fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
	
	wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");
	
	while ((wordNum==sampleNum+1)&&(!feof(fh)))
	{
		strcpy(matrix->recordInfo[recordNum].name, words[0]);
		
		for (i=0;i<sampleNum;i++)
		{
			matrix->matrix[recordNum*sampleNum+i] = atof(words[i+1]);
		}
		
		recordNum++;
		
		fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
		wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");
	}
	
	fclose(fh);
	
	free(tmpS);
	FreeWords(words, MAX_SAMPLE_NUM+1);
	
	return 1;
}

//Intersect sample ID sets in two matrics and generate new matrics with the intersected ID set. Return the number of intersected IDs
int IntersectSampleIDs(DATA_MATRIX_STRUCT *srcMatrix1, DATA_MATRIX_STRUCT *srcMatrix2, DATA_MATRIX_STRUCT *destMatrix1, DATA_MATRIX_STRUCT *destMatrix2)
{
	int i,j,k;
	int intersectSampleNum,index;
	
	destMatrix1->recordNum = srcMatrix1->recordNum;
	destMatrix2->recordNum = srcMatrix2->recordNum;
		
	intersectSampleNum = 0;
	
	for (i=0;i<srcMatrix1->sampleNum;i++)
	{
		for (j=0;j<srcMatrix2->sampleNum;j++)
		{
			if (!strcmp(srcMatrix1->sampleInfo[i].name,srcMatrix2->sampleInfo[j].name))
			{
				intersectSampleNum ++;
				break;
			}
		}
	}
	
	if ((AllocDataMatrix(destMatrix1, intersectSampleNum, srcMatrix1->recordNum)<=0)
		||(AllocDataMatrix(destMatrix2, intersectSampleNum, srcMatrix2->recordNum)<=0))
	{
		return -1;
	}
	
	memcpy(destMatrix1->recordInfo,srcMatrix1->recordInfo,srcMatrix1->recordNum*sizeof(ID_INFO_STRUCT));
	memcpy(destMatrix2->recordInfo,srcMatrix2->recordInfo,srcMatrix2->recordNum*sizeof(ID_INFO_STRUCT));
	
	index = 0;
	
	for (i=0;i<srcMatrix1->sampleNum;i++)
	{
		for (j=0;j<srcMatrix2->sampleNum;j++)
		{
			if (!strcmp(srcMatrix1->sampleInfo[i].name,srcMatrix2->sampleInfo[j].name))
			{
				memcpy(&(destMatrix1->sampleInfo[index]), &(srcMatrix1->sampleInfo[i]), sizeof(ID_INFO_STRUCT));
				memcpy(&(destMatrix2->sampleInfo[index]), &(srcMatrix2->sampleInfo[j]), sizeof(ID_INFO_STRUCT));
				
				for (k=0;k<srcMatrix1->recordNum;k++)
				{
					destMatrix1->matrix[k*intersectSampleNum+index] = srcMatrix1->matrix[k*srcMatrix1->sampleNum+i];
				}
				
				for (k=0;k<srcMatrix2->recordNum;k++)
				{
					destMatrix2->matrix[k*intersectSampleNum+index] = srcMatrix2->matrix[k*srcMatrix2->sampleNum+j];
				}
				
				index++;
			}
		}
	}
	
	return 1;
}

//Save the data matrix
int SaveDataMatrix(char *fileName, DATA_MATRIX_STRUCT *matrix)
{
	FILE *fh;
	int i,j;
	
	fh = (FILE *)fopen(fileName, "w");
	
	assert(fh!=NULL);
	
	if (!fh)
	{
		return -1;
	}
	
	fprintf(fh, "#");
	
	for (i=0;i<matrix->sampleNum;i++)
	{
		fprintf(fh, "\t%s", matrix->sampleInfo[i].name);
	}
	
	fprintf(fh,"\n");
	
	for (i=0;i<matrix->recordNum;i++)
	{
		fprintf(fh, "%s", matrix->recordInfo[i].name);
		
		for (j=0;j<matrix->sampleNum;j++)
		{
			fprintf(fh, "\t%f", matrix->matrix[i*matrix->sampleNum+j]);
		}
		
		fprintf(fh,"\n");
	}
	
	return 1;
}