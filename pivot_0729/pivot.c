//编译命令：gcc pivot.c -lm -O3 -o pivot
//运行命令：./pivot
//编译器版本：gcc 7.5.0

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>

// Calculate sum of distance while combining different pivots. Complexity : O( n^2 ) 看到这里
double SumDistance(const int k, const int n, const int dim, double* coord, int* pivots){
    double* rebuiltCoord = (double*)malloc(sizeof(double) * n * k);
    int i;
    for(i=0; i<n*k; i++){
        rebuiltCoord[i] = 0;//初始化
    }

    // Rebuild coordinates. New coordinate of one point is its distance to each pivot.
    for(i=0; i<n; i++){
        int ki;
        for(ki=0; ki<k; ki++){
            double distance = 0;
            int pivoti = pivots[ki];
            int j;
            for(j=0; j<dim; j++){
                distance += pow(coord[pivoti*dim + j] - coord[i*dim + j] ,2);// (x-y)^2
            }
            rebuiltCoord[i*k + ki] = sqrt(distance);//二维矩阵赋值
        }
    }

    // Calculate the sum of Chebyshev distance with rebuilt coordinates between every points
    double chebyshevSum = 0;
    for(i=0; i<n; i++){
        int j;
        for(j=i+1; j<n; j++){
            double chebyshev = 0;
            int ki;
            for(ki=0; ki<k; ki++){
                double dis = fabs(rebuiltCoord[i*k + ki] - rebuiltCoord[j*k + ki]);
                chebyshev = dis>chebyshev ? dis : chebyshev;
            }
            chebyshevSum += chebyshev;
        }
    }
    chebyshevSum *= 2;
    free(rebuiltCoord);

    return chebyshevSum;
}

// Recursive function Combination() : combine pivots and calculate the sum of distance while combining different pivots.
// ki  : 0 current depth of the recursion
// k   : 2 number of pivots  
// n   : 500 number of points 
// dim : 2 dimension of metric space 
// M   : 1000 number of combinations to store 
// coord  : coord coordinates of points
// pivots : temp[3] temp[0]=-1 indexes of pivots
// maxDistanceSum  : the largest M distance sum
// maxDisSumPivots : the top M pivots combinations
// minDistanceSum  : the smallest M distance sum
// minDisSumPivots : the bottom M pivots combinations
void Combination(int ki, const int k, const int n, const int dim, const int M, double* coord, int* pivots,
                 double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots){
    if(ki==k-1){
        int i;
        for(i=pivots[ki-1]+1; i<n; i++){
            pivots[ki] = i;//1->499
            struct timeval TP1, TP2;
            // Calculate sum of distance while combining different pivots.
            gettimeofday(&TP1, NULL);
            double distanceSum = SumDistance(k, n, dim, coord, pivots);
            gettimeofday(&TP2, NULL);
            if(i==n-1) printf("Using time : %f ms\n", (TP2.tv_sec - TP1.tv_sec) * 1000.0 + (TP2.tv_usec - TP1.tv_usec) / 1000.0);
            // put data at the end of array
            maxDistanceSum[M] = distanceSum;
            minDistanceSum[M] = distanceSum;
            int kj;
            for(kj=0; kj<k; kj++){
                maxDisSumPivots[M*k + kj] = pivots[kj];
            }
            for(kj=0; kj<k; kj++){
                minDisSumPivots[M*k + kj] = pivots[kj];
            }
            // sort
            int a;
            for(a=M; a>0; a--){
                if(maxDistanceSum[a] > maxDistanceSum[a-1]){
                    double temp = maxDistanceSum[a];
                    maxDistanceSum[a] = maxDistanceSum[a-1];
                    maxDistanceSum[a-1] = temp;
                    int kj;
                    for(kj=0; kj<k; kj++){
                        int temp = maxDisSumPivots[a*k + kj];
                        maxDisSumPivots[a*k + kj] = maxDisSumPivots[(a-1)*k + kj];
                        maxDisSumPivots[(a-1)*k + kj] = temp;
                    }
                }
            }
            for(a=M; a>0; a--){
                if(minDistanceSum[a] < minDistanceSum[a-1]){
                    double temp = minDistanceSum[a];
                    minDistanceSum[a] = minDistanceSum[a-1];
                    minDistanceSum[a-1] = temp;
                    int kj;
                    for(kj=0; kj<k; kj++){
                        int temp = minDisSumPivots[a*k + kj];
                        minDisSumPivots[a*k + kj] = minDisSumPivots[(a-1)*k + kj];
                        minDisSumPivots[(a-1)*k + kj] = temp;
                    }
                }
            }
        }
        return;
    }

    // Recursively call Combination() to combine pivots
    int i;
    // printf("pivots:%d\t",pivots[ki-1]+1);
    for(i=pivots[ki-1]+1; i<n; i++) {
        pivots[ki] = i;//把第一项变为0 pivots里面所含的是索引
        struct timeval TP1,TP2;
        gettimeofday(&TP1, NULL);
        Combination(ki+1, k, n, dim, M, coord, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);
        gettimeofday(&TP2, NULL);

        /** Iteration Log : pivots computed, best pivots, max distance sum, min distance sum pivots, min distance sum
        *** You can delete the logging code. **/
        if(ki==k-2){
            printf("Using time : %f ms\n", (TP2.tv_sec - TP1.tv_sec) * 1000.0 + (TP2.tv_usec - TP1.tv_usec) / 1000.0);
            int kj;
            for(kj=0; kj<k; kj++){
                printf("pivots:%d ", pivots[kj]);
            }
            for(kj=0; kj<k; kj++){
                printf("MaxPivots:%d \t", maxDisSumPivots[kj]);
            }
            printf("maxDistanceSum: %lf\t", maxDistanceSum[0]);
            for(kj=0; kj<k; kj++){
                printf("minPivots:%d \t", minDisSumPivots[kj]);
            }
            printf("minDistanceSum:%lf\n", minDistanceSum[0]);
        }
    }
}

int main(int argc, char* argv[]){
    // filename : input file namespace
    char* filename = (char*)"uniformvector-2dim-5h.txt";
    if( argc==2 ) {
        filename = argv[1];
    }  else if(argc != 1) {
        printf("Usage: ./pivot <filename>\n");
        return -1;
    }
    // M : number of combinations to store
    const int M = 1000;
    // dim : dimension of metric space
    int dim;
    // n : number of points
    int n;
    // k : number of pivots
    int k;

    // Read parameter
    FILE* file = fopen(filename, "r");
    if( file == NULL ) {
        printf("%s file not found.\n", filename);
        return -1;
    }
    fscanf(file, "%d", &dim);//2
    fscanf(file, "%d", &n);//500
    fscanf(file, "%d", &k);//2
    printf("dim = %d, n = %d, k = %d\n", dim, n, k);

    // Start timing
    struct timeval start;
    gettimeofday(&start, NULL);

    // Read Data
    double* coord = (double*)malloc(sizeof(double) * dim * n);//500*2
    int i;
    for(i=0; i<n; i++){
        int j;
        for(j=0; j<dim; j++){
            fscanf(file, "%lf", &coord[i*dim + j]);//读取所有矩阵数据
        }
    }
    fclose(file);


    // maxDistanceSum : the largest M distance sum
    // minDistanceSum : the smallest M distance sum
    double* maxDistanceSum = (double*)malloc(sizeof(double) * (M+1));//1001的空间,考虑存储空间,double 64
    double* minDistanceSum = (double*)malloc(sizeof(double) * (M + 1));
    for(i=0; i<M; i++){
        maxDistanceSum[i] = 0;
        minDistanceSum[i] = __DBL_MAX__;
    }

    // maxDisSumPivots : the top M pivots combinations
    // minDisSumPivots : the bottom M pivots combinations
    int* maxDisSumPivots = (int*)malloc(sizeof(int) * k * (M+1));//2*1001
    int* minDisSumPivots = (int*)malloc(sizeof(int) * k * (M + 1));
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k; ki++){
            maxDisSumPivots[i*k + ki] = 0;//创建了2*1001的为0的矩阵
            minDisSumPivots[i * k + ki] = 0;
        }
    }

    // temp : indexes of pivots with dummy array head
    int* temp = (int*)malloc(sizeof(int) * (k+1));//pivots矩阵合集
    temp[0] = -1;

    // Main loop. Combine different pivots with recursive function and evaluate them. Complexity : O( n^(k+2) )
    Combination(0, k, n, dim, M, coord, &temp[1], maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);

    // End timing
    struct timeval end;
    gettimeofday (&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec-start.tv_sec)*1000.0+(end.tv_usec-start.tv_usec)/1000.0);

    // Store the result
    FILE* out = fopen("result.txt", "w");
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k-1; ki++){
            fprintf(out, "%d ", maxDisSumPivots[i*k + ki]);
        }
        fprintf(out, "%d\n", maxDisSumPivots[i*k + k-1]);
    }
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k-1; ki++){
            fprintf(out, "%d ", minDisSumPivots[i*k + ki]);
        }
        fprintf(out, "%d\n", minDisSumPivots[i*k + k-1]);
    }
    fclose(file);

    // Log
    int ki;
    printf("max : ");
    for(ki=0; ki<k; ki++){
        printf("%d ", maxDisSumPivots[ki]);
    }
    printf("%lf\n", maxDistanceSum[0]);
    printf("min : ");
    for(ki=0; ki<k; ki++){
        printf("%d ", minDisSumPivots[ki]);
    }
    printf("%lf\n", minDistanceSum[0]);
    // for(i=0; i<M; i++){
        // int ki;
        // for(ki=0; ki<k; ki++){
            // printf("%d\t", maxDisSumPivots[i*k + ki]);
        // }
        // printf("%lf\n", maxDistanceSum[i]);
    // }
    // for(i=0; i<M; i++){
        // int ki;
        // for(ki=0; ki<k; ki++){
            // printf("%d\t", minDisSumPivots[i*k + ki]);
        // }
        // printf("%lf\n", minDistanceSum[i]);
    // }

    return 0;
}
