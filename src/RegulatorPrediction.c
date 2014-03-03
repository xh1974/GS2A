#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include "rngs.h"
#include "rvgs.h"
#include "words.h"
#include "math_api.h"

#define MAX_SAMPLE_NUM 1000
#define MAX_WORD_SIZE  1000
#define PERMUTATION_TIMES 1000

typedef struct
{
	char name[MAX_WORD_SIZE];
	double *values;
	int *rank;
	double *normValues;
	double score;
	double pvalue;
}GENE_DATA_STRUCT;

//Allocate the gene expression data structure
int AllocExpressionStruct(GENE_DATA_STRUCT **pExpressions, int geneNum, int sampleNum);

//Free the gene expression data structure
int FreeExpressionStruct(GENE_DATA_STRUCT **pExpressions, int geneNum);

//Read expression data from a file
int ReadExprssionFile(char *fileName, GENE_DATA_STRUCT **pExpressions, int *pGeneNum, int *pSampleNum);

//Search in gene expression data structures to match an ID. Return the index of the matched ID in the data. Called by SearchIDsFromFile()
int SearchID(char *ID, GENE_DATA_STRUCT *srcArray, int srcArrayNum);

//Search in gene expression data structures to match a list of IDs from a file. Store the pointers of matched entries in pIDMatchedData
int SearchIDsFromFile(char *fileName, GENE_DATA_STRUCT *srcArray, int srcGeneNum, GENE_DATA_STRUCT **pIDMatchedData);

//Normalize dataset using Normal Score Transformation
void NSTNormData(GENE_DATA_STRUCT *data, int geneNum, int sampleNum);

//Compute the scores and p-values for each candidate based on distance correlation
void ComputeScores(GENE_DATA_STRUCT **targetData, int targetNum, GENE_DATA_STRUCT **candData, int candNum, GENE_DATA_STRUCT *knownRegulatorData, int sampleNum, int permutationNum);

//Write to output file
int WriteToOutput(char *fileName, GENE_DATA_STRUCT **candData, int candNum);

//print command usage 
void PrintCommandUsage();

//Allocate the gene expression data structure. Called by ReadExprssionFile()
int AllocExpressionStruct(GENE_DATA_STRUCT **pExpressions, int geneNum, int sampleNum)
{
	int i;
	
	if ((geneNum<=0)||(sampleNum<=0))
	{
		return -1;
	}
	
	*pExpressions = (GENE_DATA_STRUCT *)malloc(geneNum*sizeof(GENE_DATA_STRUCT));
	
	assert(*pExpressions != NULL);
	
	for (i=0;i<geneNum;i++)
	{
		(*pExpressions)[i].values = (double *)malloc(sampleNum*sizeof(double));
		(*pExpressions)[i].rank = (int *)malloc(sampleNum*sizeof(double));
		(*pExpressions)[i].normValues = (double *)malloc(sampleNum*sizeof(double));
		
		assert(((*pExpressions)[i].values!=NULL)&&((*pExpressions)[i].rank != NULL)&&((*pExpressions)[i].normValues != NULL));
	}
	
	return geneNum;
}

//Free the gene expression data structure
int FreeExpressionStruct(GENE_DATA_STRUCT **pExpressions, int geneNum)
{
	int i;
	
	if (geneNum<=0)
	{
		return -1;
	}
	
	for (i=0;i<geneNum;i++)
	{
		free((*pExpressions)[i].values);
		free((*pExpressions)[i].rank);
		free((*pExpressions)[i].normValues);
	}
	
	free(*pExpressions);
	
	return 1;
}

//Read expression data from a file
int ReadExprssionFile(char *fileName, GENE_DATA_STRUCT **pExpressions, int *pGeneNum, int *pSampleNum)
{
	FILE *fh;
	int geneNum = 0, sampleNum = 0, wordNum = 0;
	char **words;
	char *tmpS;
	int i;
	
	fh = (FILE *)fopen(fileName, "r");
	
	assert(fh!=NULL);
	
	if (!fh)
	{
		printf("ERROR: cannot open file %s\n", fileName);
		return -1;
	}
	
	words = AllocWords(MAX_SAMPLE_NUM+1, MAX_WORD_SIZE+1);
	
	if (words == NULL)
	{
		printf("ERROR: cannot allocate memory.\n");
		return -1;
	}
	
	tmpS = (char *)malloc((MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char));
	
	assert(tmpS!=NULL);
	
	//Read the header row to get the number of samples
	fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
	
	wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");
	
	sampleNum = wordNum-1;
	
	fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
	wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");
	
	//Read through the file to get the number of genes
	while ((wordNum==sampleNum+1)&&(!feof(fh)))
	{
		geneNum++;
		
		fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
		wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");
	}
	
	fclose(fh);
	
	if (geneNum<=0)
	{
		printf("No valid data or data format error!\n");
		
		free(tmpS);
		FreeWords(words, MAX_SAMPLE_NUM+1);
		return -1;
	}
	
	//allocate memory for expression data matrix
	if (AllocExpressionStruct(pExpressions, geneNum, sampleNum)<=0)
	{
		printf("ERROR: cannot allocate memory.\n");
		free(tmpS);
		FreeWords(words, MAX_SAMPLE_NUM+1);
		return -1;
	}
	
	fh = (FILE *)fopen(fileName, "r");
	
	assert(fh!=NULL);
	
	geneNum = 0;
	
	//Read the file again to get the expression values
	
	fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
	fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
	wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");

	while ((wordNum==sampleNum+1)&&(!feof(fh)))
	{
		strcpy((*pExpressions)[geneNum].name, words[0]);
		
		for (i=0;i<sampleNum;i++)
		{
			(*pExpressions)[geneNum].values[i] = atof(words[i+1])+Uniform(-0.00000001, 0.00000001); //If two values equal, add a small random value to make them different
			(*pExpressions)[geneNum].rank[i] = -1;
			(*pExpressions)[geneNum].normValues[i] = 0.0;
		}
		
		geneNum++;
		
		fgets(tmpS, (MAX_SAMPLE_NUM+1)*(MAX_WORD_SIZE+1)*sizeof(char), fh);
		wordNum = StringToWords(words, tmpS, MAX_WORD_SIZE+1, MAX_SAMPLE_NUM+1, " \t\r\n\v\f");
	}
	
	fclose(fh);
	
	free(tmpS);
	
	*pGeneNum = geneNum;
	*pSampleNum = sampleNum;
	
	return geneNum;
}

//Search in gene expression data structures to match an ID. Return the index of the matched ID in the data. Called by SearchIDsFromFile()
int SearchID(char *ID, GENE_DATA_STRUCT *srcArray, int srcArrayNum)
{	
	int j;
	int isFound = 0;
	
	for (j=0;j<srcArrayNum;j++)
	{
		if (!strcmp(ID,srcArray[j].name))
		{
			isFound = 1;
			break;
		}
	}
	
	if (isFound)
	{
		return j;
	}
	else
	{
		return -1;
	}
}

//Search in gene expression data structures to match a list of IDs from a file. Store the pointers of matched entries in pIDMatchedData
int SearchIDsFromFile(char *fileName, GENE_DATA_STRUCT *srcArray, int srcGeneNum, GENE_DATA_STRUCT **pIDMatchedData)
{
	FILE *fh;
	char tmpS[MAX_WORD_SIZE];
	int matchedIDNum, tmpIndex;
	
	fh = (FILE *)fopen(fileName, "r");
	
	assert(fh!=NULL);
	
	if (!fh)
	{
		printf("ERROR: cannot open file %s\n", fileName);
		return -1;
	}
	
	matchedIDNum = 0;
	
	fscanf(fh,"%s",tmpS);
	
	while (!feof(fh))
	{
		tmpIndex = SearchID(tmpS, srcArray, srcGeneNum);
		
		if (tmpIndex>=0)
		{
			pIDMatchedData[matchedIDNum] = srcArray+tmpIndex;
			matchedIDNum++;
		}
		
		fscanf(fh,"%s",tmpS);
	}
	
	fclose(fh);
	
	return matchedIDNum;
}

//Normalize dataset using Normal Score Transformation
void NSTNormData(GENE_DATA_STRUCT *data, int geneNum, int sampleNum)
{
	int i;
	
	for (i=0;i<geneNum;i++)
	{
		Ranking(data[i].rank, data[i].values, sampleNum);
		NormalTransform(data[i].normValues, data[i].rank, sampleNum);
	}
}

//Compute the scores and p-values for each candidate based on distance correlation
void ComputeScores(GENE_DATA_STRUCT **targetData, int targetNum, GENE_DATA_STRUCT **candData, int candNum, GENE_DATA_STRUCT *knownRegulatorData, int sampleNum, int permutationNum)
{
	int i,j;
	double targetSignatureValues[MAX_SAMPLE_NUM];
	double regulatorValues[MAX_SAMPLE_NUM];
	double *randomScore;
	
	assert((targetNum>0)&&(candNum>0)&&(sampleNum>0)&&(permutationNum>0));
	
	//Compute the average of all targets
	
	memset(targetSignatureValues, 0, MAX_SAMPLE_NUM*sizeof(double));
	
	for (i=0;i<targetNum;i++)
	{
		for (j=0;j<sampleNum;j++)
		{
			targetSignatureValues[j] += targetData[i]->normValues[j];
		}
	}
	
	for (j=0;j<sampleNum;j++)
	{
		targetSignatureValues[j] /= targetNum;
	}
	
	//Compute distance correlation for each candidate
	memcpy(regulatorValues, knownRegulatorData->normValues, sampleNum*sizeof(double));
	
	printf("Processing candidates.....\n");
	
	for (i=0;i<candNum;i++)
	{
		memcpy(regulatorValues+sampleNum, candData[i]->normValues, sampleNum*sizeof(double));
		candData[i]->score = ComputeDistanceCorrelation(regulatorValues, targetSignatureValues, 2, sampleNum);
		
		if (i%(candNum/100)==0)
		{
			printf("%d percent of candidates processed. \r", i/(candNum/100));
		}
	}
	
	//Generate random score based on permutation
	
	printf("Permutation ......\n");
	
	randomScore = (double *)malloc(permutationNum*sizeof(double));
	
	assert(randomScore!=NULL);
	
	for (i=0;i<permutationNum;i++)
	{
		PermuteFloatArrays(regulatorValues+sampleNum, sampleNum);
		
		randomScore[i] = ComputeDistanceCorrelation(regulatorValues, targetSignatureValues, 2, sampleNum);
		
		if (i%(permutationNum/100)==0)
		{
			printf("%d percent permutation passes finished. \r", i/(permutationNum/100));
		}
	}
	
	//Compute pvalues for each candidate
	
	QuicksortF(randomScore, 0, permutationNum-1);
	
	for (i=0;i<candNum;i++)
	{
		candData[i]->pvalue = (permutationNum-0.5-bTreeSearchingF(candData[i]->score, randomScore, 0, permutationNum-1))/permutationNum;
	}
	
	free(randomScore);
}

//Write to output file
int WriteToOutput(char *fileName, GENE_DATA_STRUCT **candData, int candNum)
{
	FILE *fh;
	int i;
	
	fh = (FILE *)fopen(fileName, "w");
	
	if (!fh)
	{
		return 0;
	}
	
	for (i=0;i<candNum;i++)
	{
		fprintf(fh, "%s\t%f\t%f\n", candData[i]->name, candData[i]->score, candData[i]->pvalue);
	}
	
	fclose(fh);
	
	return 1;
}

//print command usage 
void PrintCommandUsage()
{
	//print the options of the command
	printf("RegulatorMiner - Predicting upstream regulators or co-regulators of a known regulator based on distance correlation.\n");
	printf("usage:\n");
	printf("-d <expression data file>\n");
	printf("-t <target gene id file>\n");
	printf("-c <candidate gene id file>\n");
	printf("-r <name of known regulator>\n");
	printf("-o <output file>\n");
	printf("example:\n");
	printf("RegulatorMiner -d BRCA_sample_expression.txt -t ER_target_ID.txt -c candidate_ID.txt -r ESR1 -o output.txt \n");
}

int main (int argc, const char * argv[]) 
{
	char expressionFileName[1000], targetIDFileName[1000], candidateIDFileName[1000], knownRegulatorID[1000], outputFileName[1000];
	GENE_DATA_STRUCT *expressions;
	GENE_DATA_STRUCT **targetData;
	GENE_DATA_STRUCT **candData;
	GENE_DATA_STRUCT *knownRegulatorData;
	int allGeneNum, targetGeneNum, candGeneNum, sampleNum, tmpIndex;
	int i;
	
	
	//Parse the command line
	if (argc == 1)
	{
		PrintCommandUsage();
		return -1;
	}
	
	expressionFileName[0] = 0;
	targetIDFileName[0] = 0;
	candidateIDFileName[0] = 0;
	knownRegulatorID[0] = 0;
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
			strcpy(candidateIDFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-r")==0)
		{
			strcpy(knownRegulatorID, argv[i]);
		}
		if (strcmp(argv[i-1], "-o")==0)
		{
			strcpy(outputFileName, argv[i]);
		}
	}
	
	if ((expressionFileName[0]==0)||(targetIDFileName[0]==0)||(candidateIDFileName[0]==0)||(knownRegulatorID[0]==0)||(outputFileName[0]==0))
	{
		printf("Command error!\n");
		PrintCommandUsage();
		return -1;
	}
	
	PlantSeeds(123456);
	
	//read expression data
	if (ReadExprssionFile(expressionFileName, &expressions, &allGeneNum, &sampleNum)<=0)
	{
		printf("ERROR: cannot open %s or incorrect format!\n", expressionFileName);
		return -1;
	}
	else
	{
		printf("%d genes in the expression data.\n", allGeneNum);
	}
	
	tmpIndex = SearchID(knownRegulatorID, expressions, allGeneNum);
	
	if (tmpIndex<0)
	{
		printf("ERROR: cannot find the known regulator ID in expression data!\n");
		FreeExpressionStruct(&expressions, allGeneNum);
		return -1;
	}
	
	knownRegulatorData = expressions + tmpIndex;
					  
	targetData = (GENE_DATA_STRUCT **)malloc(allGeneNum*sizeof(GENE_DATA_STRUCT *));
	candData = (GENE_DATA_STRUCT **)malloc(allGeneNum*sizeof(GENE_DATA_STRUCT *));
	
	assert((targetData!=NULL)&&(candData!=NULL));
	
	if ((targetData==NULL)||(candData==NULL))
	{
		printf("ERROR: cannot allocate memory!\n");
		FreeExpressionStruct(&expressions, allGeneNum);
		return -1;
	}
	
	//Get target gene data
	targetGeneNum = SearchIDsFromFile(targetIDFileName, expressions, allGeneNum, targetData);
	
	if (targetGeneNum<=0)
	{
		printf("no target ID found in expression data!\n");
		FreeExpressionStruct(&expressions, allGeneNum);
		free(targetData);
		free(candData);
		return -1;
	}
	else
	{
		printf("%d target genes found in expression data.\n", targetGeneNum);
	}
	
	//Get candidate gene data
	candGeneNum = SearchIDsFromFile(candidateIDFileName, expressions, allGeneNum, candData);
	
	if (candGeneNum<=0)
	{
		printf("no candidate ID found in expression data!\n");
		FreeExpressionStruct(&expressions, allGeneNum);
		free(targetData);
		free(candData);
		return -1;
	}
	else
	{
		printf("%d candidate genes found in expression data.\n", candGeneNum);
	}
	
	NSTNormData(expressions, allGeneNum, sampleNum);
	
	ComputeScores(targetData, targetGeneNum, candData, candGeneNum, knownRegulatorData, sampleNum, PERMUTATION_TIMES);
	
	if (!WriteToOutput(outputFileName, candData, candGeneNum))
	{
		printf("ERROR: cannot open %s!\n", outputFileName);
		FreeExpressionStruct(&expressions, allGeneNum);
		free(targetData);
		free(candData);
		return -1;
	}
	
	FreeExpressionStruct(&expressions, allGeneNum);
	free(targetData);
	free(candData);
	
	return 0;
}
