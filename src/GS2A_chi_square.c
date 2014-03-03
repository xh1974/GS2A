#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <float.h>
#include "rngs.h"
#include "rvgs.h"
#include "math_api.h"
#include "dataMatrix.h"

#define PERMUTATION_NUM 1000

typedef struct
{
	ID_INFO_STRUCT *id;
	double score;
	double pValue;
}CANDIDATE_SCORE_STRUCT;

//Search in gene expression data structures to mark a list of IDs in a file.
int MarkIDs(char *fileName, DATA_MATRIX_STRUCT *data);

//Compute GS2A score for a candidate feature
double ComputeGS2AScore(DATA_MATRIX_STRUCT *expressionMatrix, 
						double *feature, 
						int featureSize,
						char *maskedID);

//Compute scores for all candidates and store the values in candidate score structure
int ComputeScoreMain(DATA_MATRIX_STRUCT *expressionMatrix, DATA_MATRIX_STRUCT *candidateMatrix, CANDIDATE_SCORE_STRUCT *candidateScores);

//Compute p-values for candidates based on permutation
int ComputePermutationP(DATA_MATRIX_STRUCT *expressionMatrix, DATA_MATRIX_STRUCT *candidateMatrix, CANDIDATE_SCORE_STRUCT *candidateScores, int candidateNum, int permutationNum);

//Write to output file
int WriteToOutput(char *fileName, CANDIDATE_SCORE_STRUCT *candidateScores, int candNum);

//print command usage 
void PrintCommandUsage();

//Search in gene expression data structures to mark a list of IDs in a file. Return number of matched ID
int MarkIDs(char *fileName, DATA_MATRIX_STRUCT *data)
{
	FILE *fh;
	char tmpS[MAX_WORD_SIZE];
	int matchedIDNum;
	int i;
	
	fh = (FILE *)fopen(fileName, "r");
	
	assert(fh!=NULL);
	
	if (!fh)
	{
		printf("ERROR: cannot open file %s\n", fileName);
		return -1;
	}
	
	for (i=0;i<data->recordNum;i++)
	{
		data->recordInfo[i].flag = 0;
	}
	
	matchedIDNum = 0;
	
	fscanf(fh,"%s",tmpS);
	
	while (!feof(fh))
	{
		for (i=0;i<data->recordNum;i++)
		{
			if (!strcmp(data->recordInfo[i].name,tmpS))
			{
				data->recordInfo[i].flag = 1;
				matchedIDNum++;
				break;
			} 
		}
		
		fscanf(fh,"%s",tmpS);
	}
	
	fclose(fh);
	
	return matchedIDNum;
}

//Compute GS2A score for a candidate feature
double ComputeGS2AScore(DATA_MATRIX_STRUCT *expressionMatrix, 
					 double *feature, 
					 int featureSize,
					 char *maskedID)
{
	int i;
	int geneNum = expressionMatrix->recordNum;
	double *targetValues, *nonTargetValues;
	int targetNum = 0, nonTargetNum = 0;
	double targetSquareSum, nonTargetMean, nonTargetStdev;
		
	assert(featureSize==expressionMatrix->sampleNum);
	assert(featureSize>1);
	assert(geneNum>1);
	
	if (!((featureSize==expressionMatrix->sampleNum)&&(featureSize>1)&&(geneNum>1)))
	{
		return 1;
	}
	
	targetValues = (double *)malloc(geneNum*sizeof(double));
	nonTargetValues = (double *)malloc(geneNum*sizeof(double));
	
	assert((targetValues!=NULL)&&(nonTargetValues!=NULL));
	
	//Compute Correlation values
	for (i=0;i<geneNum;i++)
	{
		if (!strcmp(expressionMatrix->recordInfo[i].name, maskedID))
		{
			continue;
		}
		
		if (expressionMatrix->recordInfo[i].flag)
		{
			targetValues[targetNum] = PearsonCorrel(feature, &(expressionMatrix->matrix[i*featureSize]), featureSize);
			targetNum++;
		}
		else
		{
			nonTargetValues[nonTargetNum] = PearsonCorrel(feature, &(expressionMatrix->matrix[i*featureSize]), featureSize);
			nonTargetNum++;
		}
	}
	
	assert((targetNum>0)&&(nonTargetNum>0));
	
	if (!((targetNum>0)&&(nonTargetNum>0)))
	{
		free(targetValues);
		free(nonTargetValues);
		return 0;
	}
	
	nonTargetMean = 0;
	
	for (i=0;i<nonTargetNum;i++)
	{
		nonTargetMean += nonTargetValues[i];
	}
	
	nonTargetMean /= nonTargetNum;
	
	nonTargetStdev = 0;
	
	for (i=0;i<nonTargetNum;i++)
	{
		nonTargetStdev += (nonTargetValues[i]-nonTargetMean)*(nonTargetValues[i]-nonTargetMean);
	}
	
	nonTargetStdev = sqrt(nonTargetStdev/nonTargetNum);

	targetSquareSum = 0;

	for (i=0;i<targetNum;i++)
	{
		targetSquareSum += (targetValues[i]-nonTargetMean)*(targetValues[i]-nonTargetMean)/nonTargetStdev/nonTargetStdev;
	}
	
	free(targetValues);
	free(nonTargetValues);
	
	return targetSquareSum/(targetNum-1)-0.5;
}
	
//Compute p-values for candidates based on permutation
int ComputePermutationP(DATA_MATRIX_STRUCT *expressionMatrix, DATA_MATRIX_STRUCT *candidateMatrix, CANDIDATE_SCORE_STRUCT *candidateScores, int candidateNum, int permutationNum)
{
	int i;	
	double *randScore;
	int sampleNum = candidateMatrix->sampleNum;
	double *tmpFeature;
	int tmpIndex;
	
	assert(expressionMatrix->sampleNum==candidateMatrix->sampleNum);
	
	if (expressionMatrix->sampleNum!=candidateMatrix->sampleNum)
	{
		return -1;
	}
	
	randScore = (double *)malloc(permutationNum*sizeof(double));
	tmpFeature = (double *)malloc(sampleNum*sizeof(double));
	
	assert((randScore!=NULL)&&(tmpFeature!=NULL));
	
	for (i=0;i<permutationNum;i++)
	{
		tmpIndex = (int)Uniform(0, candidateMatrix->recordNum);
		tmpIndex = tmpIndex<0?0:(tmpIndex>=candidateMatrix->recordNum?candidateMatrix->recordNum-1:tmpIndex);
		
		memcpy(tmpFeature, candidateMatrix->matrix+tmpIndex*sampleNum, sampleNum*sizeof(double));
		PermuteFloatArrays(tmpFeature,sampleNum);
		
		randScore[i] = fabs(ComputeGS2AScore(expressionMatrix, tmpFeature, sampleNum, ""));
	}
	
	QuicksortF(randScore, 0, permutationNum-1);
	
	for (i=0;i<candidateNum;i++)
	{
		candidateScores[i].pValue = (double)(permutationNum-1-bTreeSearchingF(fabs(candidateScores[i].score), randScore, 0, permutationNum-1))/permutationNum;
	}

	free(randScore);
	free(tmpFeature);
	return 1;
}

//Compute scores for all candidates and store the values in candidate score structure
int ComputeScoreMain(DATA_MATRIX_STRUCT *expressionMatrix, DATA_MATRIX_STRUCT *candidateMatrix, CANDIDATE_SCORE_STRUCT *candidateScores)
{
	int i;
	double *weights;

	assert(expressionMatrix->sampleNum==candidateMatrix->sampleNum);
	
	if (expressionMatrix->sampleNum!=candidateMatrix->sampleNum)
	{
		return -1;
	}
	
	weights = (double *)malloc(expressionMatrix->recordNum*sizeof(double));
	
	assert(weights!=NULL);
	
	for (i=0;i<expressionMatrix->recordNum;i++)
	{
		weights[i] = 1.0;
	}
	
	for (i=0;i<candidateMatrix->recordNum;i++)
	{
		candidateScores[i].score = ComputeGS2AScore(expressionMatrix, 
						 candidateMatrix->matrix+i*(candidateMatrix->sampleNum), 
						 candidateMatrix->sampleNum, 
						 candidateMatrix->recordInfo[i].name);
		candidateScores[i].id = candidateMatrix->recordInfo+i;
		candidateScores[i].pValue = 1;
	}
	
	free(weights);
	
	return 1;
}

//Write to output file
int WriteToOutput(char *fileName, CANDIDATE_SCORE_STRUCT *candidateScores, int candNum)
{
	FILE *fh;
	int i;
	
	fh = (FILE *)fopen(fileName, "w");
	
	if (!fh)
	{
		return 0;
	}
	
	fprintf(fh, "ID\tisSignature\tscore\tpValue\n");
	
	for (i=0;i<candNum;i++)
	{
		fprintf(fh, "%s\t%d\t%f\t%f\n", 
				candidateScores[i].id->name,
				candidateScores[i].id->flag,
				candidateScores[i].score,
				candidateScores[i].pValue);
	}
	
	fclose(fh);
	
	return 1;
}

//print command usage 
void PrintCommandUsage()
{
	//print the options of the command
	printf("GS2A - Gene Signature Association Analysis.\n");
	printf("usage:\n");
	printf("-d <expression data file>\n");
	printf("-t <target gene id file>\n");
	printf("-c <candidate data file>\n");
	printf("-o <output file>\n");
	printf("example:\n");
	printf("GS2A -d expression.txt -t target.txt -c candidate.txt -o output.txt \n");
}

int main (int argc, const char * argv[]) 
{
	char expressionFileName[1000], targetIDFileName[1000], candidateFileName[1000], outputFileName[1000];
	DATA_MATRIX_STRUCT expressions;
	DATA_MATRIX_STRUCT candidate;
	DATA_MATRIX_STRUCT expressionTrimmed;
	DATA_MATRIX_STRUCT candidateTrimmed;
	CANDIDATE_SCORE_STRUCT *candScores;
	int matchedIDNum;
	int i;
	
	//Parse the command line
	if (argc == 1)
	{
		PrintCommandUsage();
		return -1;
	}
	
	expressionFileName[0] = 0;
	targetIDFileName[0] = 0;
	candidateFileName[0] = 0;
	outputFileName[0] = 0;
	
	for (i=2;i<argc;i++)
	{
		if (strcmp(argv[i-1], "-d")==0)
		{
			strcpy(expressionFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-t")==0)
		{
			strcpy(targetIDFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-c")==0)
		{
			strcpy(candidateFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-o")==0)
		{
			strcpy(outputFileName, argv[i]);
		}
	}
	
	if ((expressionFileName[0]==0)||(targetIDFileName[0]==0)||(candidateFileName[0]==0)||(outputFileName[0]==0))
	{
		printf("Command error!\n");
		PrintCommandUsage();
		return -1;
	}
	
	PlantSeeds(123456);
	
	//Read expression data
	
	if (ReadDataMatrix(expressionFileName, &expressions)<=0)
	{
		printf("ERROR: cannot open %s or incorrect format!\n", expressionFileName);
		return -1;
	}
	else
	{
		printf("%d genes and %d samples in expression data\n", expressions.recordNum, expressions.sampleNum);
	}
	
	//read candidate data
	
	if (ReadDataMatrix(candidateFileName, &candidate)<=0)
	{
		printf("ERROR: cannot open %s or incorrect format!\n", candidateFileName);
		FreeDataMatrix(&expressions);
		return -1;
	}
	else
	{
		printf("%d records and %d samples in candidate data\n", candidate.recordNum, candidate.sampleNum);
	}
	
	//read ID data
	
	matchedIDNum = MarkIDs(targetIDFileName, &expressions);
	
	printf("%d genes in expression data are signitures\n", matchedIDNum);
	
	if (matchedIDNum <=0)
	{
		printf("no signature gene found in expression data!\n");
		FreeDataMatrix(&expressions);
		FreeDataMatrix(&candidate);
		
		return -1;
	}
	
	MarkIDs(targetIDFileName, &candidate);
	
	//intersect expression data and candidate data by samples
	
	if (IntersectSampleIDs(&expressions, &candidate, &expressionTrimmed, &candidateTrimmed)<=0)
	{
		printf("Failed in matching samples between expression data and candidate data.");
		FreeDataMatrix(&expressions);
		FreeDataMatrix(&candidate);
		
		return -1;
	}
	else
	{
		printf("%d samples in the intersaction of expression dataset and candidate dataset.\n", candidateTrimmed.sampleNum);
	}
	
	candScores = (CANDIDATE_SCORE_STRUCT *)malloc((candidateTrimmed.recordNum)*sizeof(CANDIDATE_SCORE_STRUCT));
	
	assert(candScores!=NULL);
	
	printf("Computing GS2A scores......\n");
	
	ComputeScoreMain(&expressionTrimmed, &candidateTrimmed, candScores);
	
	printf("Permutation......\n");
	
	ComputePermutationP(&expressionTrimmed, &candidateTrimmed, candScores, candidateTrimmed.recordNum, PERMUTATION_NUM);
	
	if (!WriteToOutput(outputFileName, candScores, candidateTrimmed.recordNum))
	{
		printf("Cannot write to %s!\n", outputFileName);
	}
	
	FreeDataMatrix(&expressions);
	FreeDataMatrix(&candidate);
	FreeDataMatrix(&expressionTrimmed);
	FreeDataMatrix(&candidateTrimmed);
	free(candScores);
	
	printf("Finished.\n");
	
	return 0;
}
