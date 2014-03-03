/*
 *  dataMatrixIO.h
 *  processing data matrix
 *
 *  Created by Han Xu on 4/9/13.
 *  Copyright 2013 Dana Farber Cancer Institute. All rights reserved.
 *
 */

#define MAX_SAMPLE_NUM 10000
#define MAX_WORD_SIZE  255

typedef struct
{
	char name[MAX_WORD_SIZE];
	int flag;
}ID_INFO_STRUCT;

typedef struct
{
	ID_INFO_STRUCT *sampleInfo;
	ID_INFO_STRUCT *recordInfo;
	double *matrix;
	int sampleNum;
	int recordNum;
}DATA_MATRIX_STRUCT;

//allocate memory for data matrix
int AllocDataMatrix(DATA_MATRIX_STRUCT *matrix, int sampleNum, int recordNum);

//Free memory for data matrix
void FreeDataMatrix(DATA_MATRIX_STRUCT *matrix);

//Read the data matrix from a file
int ReadDataMatrix(char *fileName, DATA_MATRIX_STRUCT *matrix);

//Intersect sample ID sets in two matrics and generate new matrics with the intersected ID set. Return the number of intersected IDs
int IntersectSampleIDs(DATA_MATRIX_STRUCT *srcMatrix1, DATA_MATRIX_STRUCT *srcMatrix2, DATA_MATRIX_STRUCT *destMatrix1, DATA_MATRIX_STRUCT *destMatrix2);

//Save the data matrix
int SaveDataMatrix(char *fileName, DATA_MATRIX_STRUCT *matrix);
