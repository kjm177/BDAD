#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/*** Skeleton for Lab 1 ***/

/***** Globals ******/
float **a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float err; /* The absolute relative error */
int num = 0;  /* number of unknowns */


/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */

/********************************/



/* Function definitions: functions are ordered alphabetically ****/
/*****************************************************************/

/* 
   Conditions for convergence (diagonal dominance):
   1. diagonal element >= sum of all other elements of the row
   2. At least one diagonal element > sum of all other elements of the row
 */
void check_matrix()
{
  int bigger = 0; /* Set to 1 if at least one diag element > sum  */
  int i, j;
  float sum = 0;
  float aii = 0;
  
  for(i = 0; i < num; i++)
  {
    sum = 0;
    aii = fabs(a[i][i]);
    
    for(j = 0; j < num; j++)
       if( j != i)
	 sum += fabs(a[i][j]);
       
    if( aii < sum)
    {
      printf("The matrix will not converge.\n");
      exit(1);
    }
    
    if(aii > sum)
      bigger++;
    
  }
  
  if( !bigger )
  {
     printf("The matrix will not converge\n");
     exit(1);
  }
}


float* allocateFloatArray(int n)
{
    float* temp = (float *)malloc(n * sizeof(float)); 
    if(!temp)
  	{
		printf("ERROR! Cannot allocate memory\n");
		exit(1);
  	}
	return temp;
}


/******************************************************/
/* Read input from file */
/* After this function returns:
 * a[][] will be filled with coefficients and you can access them using a[i][j] for element (i,j)
 * x[] will contain the initial values of x
 * b[] will contain the constants (i.e. the right-hand-side of the equations
 * num will have number of variables
 * err will have the absolute error that you need to reach
 */
 
 
 
void get_input(char filename[])
{
  FILE * fp;
  int i,j;  
 
  fp = fopen(filename, "r");
  if(!fp)
  {
    printf("Cannot open file %s\n", filename);
    exit(1);
  }

 fscanf(fp,"%d ",&num);
 fscanf(fp,"%f ",&err);

 /* Now, time to allocate the matrices and vectors */
 a = (float**)malloc(num * sizeof(float*));
 if( !a)
  {
	printf("Cannot allocate a!\n");
	exit(1);
  }

	for(i = 0; i < num; i++) 
		a[i] = allocateFloatArray(num);
 
	x = allocateFloatArray(num);
	b = allocateFloatArray(num);

 /* Now .. Filling the blanks */ 

 /* The initial values of Xs */
 for(i = 0; i < num; i++)
	fscanf(fp,"%f ", &x[i]);
 
 for(i = 0; i < num; i++)
 {
   for(j = 0; j < num; j++)
     fscanf(fp,"%f ",&a[i][j]);
   
   /* reading the b element */
   fscanf(fp,"%f ",&b[i]);
 }
 
 fclose(fp); 

}

int checkError(float *current, int myRank)
{
    if(fabs(current[myRank] - x[myRank])/current[myRank] > err)
            return 0;
	return 1;
}


int validError(float* current)
{
	int i;
    for(i = 0 ; i < num ; i++)
    {
        if(checkError(current, i) == 0)
            return 0;
    }
    return 1;
}


void solveEquation(float* y, int myRank)
{
    float sum = 0;
	int i;
    for(i = 0 ; i < num ; i++)
    {
        if(i != myRank)
        {
            sum += a[myRank][i] * x[i];
        }
    }
    y[myRank] = (b[myRank] - sum)/a[myRank][myRank];
}


/************************************************************/


int main(int argc, char *argv[])
{

	int nit = 0; /* number of iterations */
	FILE * fp;
	char output[100] ="";
	
	int i, j = 0;
	int myRank, numOfProcesses, numOfDivision, currentError;
	

	if( argc != 2)
	{
	printf("Usage: ./gsref filename\n");
	exit(1);
	}

	/* Read the input file and fill the global data structure above */ 
	get_input(argv[1]);

	/* Check for convergence condition */
	/* This function will exit the program if the coffeicient will never converge to 
	* the needed absolute error. 
	* This is not expected to happen for this programming assignment.
	*/
	check_matrix();

 
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Status stat;
 
	if(numOfProcesses == 1)
	{
		while(1)
		{
			float *y;
			y = allocateFloatArray(num);
			for(i=0;i<num;i++)
				solveEquation(y,i);
			nit++;
			
			if(validError(y) != 0)
			{
				for(i = 0; i < num; i++)
			       	x[i] = y[i];	
				break;
			}
			for(i = 0; i < num; i++)
				x[i] = y[i];	
			
			free(y);
		}
	}
	else
	{
		numOfDivision = ceil((double)num/(double)(numOfProcesses-1));
		if(myRank == 0)
		{
			while(j < num)
			{
				j = 0;
				float *tempX = (float*)malloc(sizeof(float) * (num+1));
				memcpy(tempX, x, sizeof(float) * num);
				tempX[num] = 0;
				MPI_Bcast((void *)tempX, num+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
			
				if(numOfProcesses > num)
				{
					printf("CAUTION! Number of processes > number of unknowns. \n This case is not part of lab 2.\n");
					exit(1);
				}
				else
				{
					j = 0;
					for(i = 0; i < numOfProcesses-1; i++)
					{
						float *updatedX = allocateFloatArray(num);
						MPI_Recv((void *)updatedX, num, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&stat);
						int k;
						for(k = (stat.MPI_SOURCE -1)*numOfDivision; k < (stat.MPI_SOURCE -1) * numOfDivision + numOfDivision ; k++)
						{
								if(k >= num)
									break;
								x[k] = updatedX[k];
								if(stat.MPI_TAG == 1)
									j++;
						}
					}
					nit++;
				}
			}
			float *tempX = (float*)malloc(sizeof(float) * (num+1));
			tempX[num] = 1;
			MPI_Bcast((void *)tempX, num+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		}
		else
		{
			while(1)
			{
				float* tempX = (float*)malloc(sizeof(float) * (num+1));
				MPI_Bcast((void *)tempX, num+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
				if (tempX[num] == 1)
					break;
				memcpy(x, tempX, sizeof(float)*num);
				
				if(numOfProcesses > num)
				{
					printf("CAUTION! Number of processes > number of unknowns. \n This case is not part of lab 2.\n");
					exit(1);
				}
				else if(numOfProcesses <= num)
				{
					numOfDivision = ceil((double)num/(double)(numOfProcesses-1));
					int k;
					int tempError = 1;
					float *z = allocateFloatArray(num);
					for(k = (myRank - 1) * numOfDivision ; k < (myRank - 1) * numOfDivision + numOfDivision ; k++)
					{
						if(k >= num)
							break;
						float *y = allocateFloatArray(num);
						solveEquation(y, k);
						
						if(checkError(y, k) == 0)
							tempError = 0;
						
						z[k] = y[k];						
					}
					MPI_Send((void *)z, num, MPI_FLOAT, 0, tempError, MPI_COMM_WORLD);
				}	
			}
		}	
	}

 /* Writing results to file */
	if(myRank == 0)
	{
		sprintf(output,"%d.sol",num);
		fp = fopen(output,"w");
		if(!fp)
		{
			printf("Cannot create the file %s\n", output);
			exit(1);
		}

		for( i = 0; i < num; i++)
			fprintf(fp,"%f\n",x[i]);

		printf("total number of iterations: %d\n", nit);

		fclose(fp);
	}

	MPI_Finalize();
	exit(0);

}