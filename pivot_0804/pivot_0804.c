/*******
    run.sh:
    #!/bin/bash
    #SBATCH --job-name=Pivot0804
    #SBATCH --partition=IPCC //如果在学校超算为 cu
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=32
    #SBATCH --error=%j.err
    #SBATCH --output=%j.out
    icpc pivot_0804.c -Ofast -lm -qopenmp -o pivot_0804
    ./pivot_0804

    Add batch work: sbatch run.sh
*******/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<omp.h>

void SumP2PDistance(double* P2PDist, double* coord, int n, int dim) {
    int i;
#pragma omp parallel for num_threads(32)
    for (i = 0; i < n; i++) {
        int ki;
        for (ki = 0; ki < n; ki++) {
            double distance = 0;
            int j;
            for (j = 0; j < dim; j++) { 
                distance += pow(coord[ki * dim + j] - coord[i * dim + j], 2);
            }
            //Array Content: Distance[i,ki]
            P2PDist[i * n + ki] = sqrt(distance);
            //Testlog:
            //if(i==0) printf("point1:%d\tpoint2:%d\tdistance:%lf\n", i, ki, P2PDist[i * (n - 1) + ki - subnum]);
        }
    }
    return;
}

double SumDistance(const int k, const int n, const int dim, double* coord, int* pivots, double* P2PDist) {
    int i;
    double chebyshevSum = 0;
#pragma omp parallel for reduction(+:chebyshevSum) num_threads(32)
    for (i = 0; i < n; i++) {
        int j;
        for (j = i + 1; j < n; j++) {   //openmp reduction sum
            double chebyshev = 0;
            int ki;
            for (ki = 0; ki < k; ki++) {    //  MAX(  |d(xi,pki)-d(xj,pki)| )   openmp reduction max
                //double dis = fabs(getP2PDistance(i,pivots[ki],P2PDist,n) - getP2PDistance(j, pivots[ki], P2PDist, n));
                double dis = fabs(P2PDist[i * n + pivots[ki]] - P2PDist[j * n + pivots[ki]]);
                chebyshev = dis > chebyshev ? dis : chebyshev;
            }
            chebyshevSum += chebyshev;
        }
    }
    chebyshevSum *= 2;
    return chebyshevSum;
}

//Calculate sum of distance while combining different pivots. Complexity : O( n^2 )
//double SumDistance(const int k, const int n, const int dim, double* coord, int* pivots, double* p2pdist) {
//    double* rebuiltcoord = (double*)malloc(sizeof(double) * n * k);
//    int i;
//    for (i = 0; i < n * k; i++) {
//        rebuiltcoord[i] = 0;//初始化
//    }
//
//    // rebuild coordinates. new coordinate of one point is its distance to each pivot.
//    for (i = 0; i < n; i++) {
//        int ki;
//        for (ki = 0; ki < k; ki++) {
//            int pivoti = pivots[ki];
//            int j;
//            double distance = 0;
//            for (j = 0; j < dim; j++) {
//                distance += pow(coord[pivoti * dim + j] - coord[i * dim + j], 2);// (pivoti-i)^2
//            }
//            rebuiltcoord[i * k + ki] = sqrt(distance);//二维矩阵赋值
//
//            //if (pivoti==0) printf("point1:%d\tpoint2:%d\tdistance:%lf\n", i, pivoti, distance);
//            //if (i > pivoti) {
//            ////   rebuiltcoord[i * k + ki] = p2pdist[i * (n - 1) + pivoti - sum(i)];
//            //}else if(i < pivoti) {
//            //    rebuiltcoord[i * k + ki] = p2pdist[pivoti * (n - 1) + i - sum(pivoti)];
//            //}
//            //else {
//            //    rebuiltcoord[i * k + ki] = 0;
//            //}
//            //rebuiltcoord[i * k + ki] = p2pdist[i * (n - 1) + pivoti - subnum];
//            //teslog
//            //if (i == 0) printf("point1:%d\tpoint2:%d\tdistance:%lf\n", i, ki, distance);
//
//            //下面可以直接使用
//            /*int subsum = sum(min(i, pivoti));
//            if (i == pivoti) rebuiltcoord[i * k + ki] = 0;
//            else rebuiltcoord[i * k + ki] = p2pdist[min(i, pivoti) * (n - 1) + max(i, pivoti) - subsum];*/
//        }
//    }
//    //gettimeofday(&tp2, null);
//    //printf("rebuild using time : %f ms\n", (tp2.tv_sec - tp1.tv_sec) * 1000.0 + (tp2.tv_usec - tp1.tv_usec) / 1000.0);
//   ////for (i = 0; i < k; i++) {
//   ////    int pivoti = pivots[i];
//   ////    int ki;
//   ////    int subnum = sum(pivoti);
//   ////    //ki<pivoti
//   ////    for (ki = 0; ki < pivoti; ki++) {
//   ////        rebuiltcoord[ki * k + ki] = p2pdist[pivoti * (n - 1) + i - sum(pivoti)];
//   ////    }
//   ////    //ki=pivoti
//   ////    rebuiltcoord[i * k + ki] = 0;
//   ////    //ki>pivoti
//   ////    for (ki = pivoti + 1; ki < n; ki++) {
//   ////        rebuiltcoord[i * k + ki] = p2pdist[i * (n - 1) + pivoti - sum(i)];
//   ////    }
//   ////}
//
//   // calculate the sum of chebyshev distance with rebuilt coordinates between every points
//    double chebyshevsum = 0;
//    #pragma omp parallel for reduction(+:chebyshevsum) num_threads(5)
//    for (i = 0; i < n; i++) {
//        int j;
//        for (j = i + 1; j < n; j++) {
//            double chebyshev = 0;
//            int ki;
//            for (ki = 0; ki < k; ki++) {
//                double dis = fabs(rebuiltcoord[i * k + ki] - rebuiltcoord[j * k + ki]); //可以优化整合？
//                chebyshev = dis > chebyshev ? dis : chebyshev;
//            }
//            chebyshevsum += chebyshev;
//        }
//    }
//    chebyshevsum *= 2;
//    free(rebuiltcoord);
//
//    return chebyshevsum;
//}

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
    double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots, double* P2PDist) {
    //end of recur
    if (ki == k - 1) {
        int i;
        for (i = pivots[ki - 1] + 1; i < n; i++) {
            pivots[ki] = i;//1->499
            //struct timeval TP1, TP2;
            // Calculate sum of distance while combining different pivots.
            //gettimeofday(&TP1, NULL);
            double distanceSum = SumDistance(k, n, dim, coord, pivots, P2PDist);
            //gettimeofday(&TP2, NULL);
            //if (i == n - 1) printf("SumDistance Using time : %f ms\n", (TP2.tv_sec - TP1.tv_sec) * 1000.0 + (TP2.tv_usec - TP1.tv_usec) / 1000.0);
            // put data at the end of array
            maxDistanceSum[M] = distanceSum;
            minDistanceSum[M] = distanceSum;
            int kj;
            for (kj = 0; kj < k; kj++) {
                maxDisSumPivots[M * k + kj] = pivots[kj];
                minDisSumPivots[M * k + kj] = pivots[kj];
            }
            // sort
            gettimeofday(&TP1, NULL);
            int a;
            for (a = M; a > 0; a--) {
                if (maxDistanceSum[a] > maxDistanceSum[a - 1]) {
                    double temp = maxDistanceSum[a];
                    maxDistanceSum[a] = maxDistanceSum[a - 1];
                    maxDistanceSum[a - 1] = temp;
                    int kj;
                    for (kj = 0; kj < k; kj++) {
                        int temp = maxDisSumPivots[a * k + kj];
                        maxDisSumPivots[a * k + kj] = maxDisSumPivots[(a - 1) * k + kj];
                        maxDisSumPivots[(a - 1) * k + kj] = temp;
                    }
                }
                if (minDistanceSum[a] < minDistanceSum[a - 1]) {
                    double temp = minDistanceSum[a];
                    minDistanceSum[a] = minDistanceSum[a - 1];
                    minDistanceSum[a - 1] = temp;
                    int kj;
                    for (kj = 0; kj < k; kj++) {
                        int temp = minDisSumPivots[a * k + kj];
                        minDisSumPivots[a * k + kj] = minDisSumPivots[(a - 1) * k + kj];
                        minDisSumPivots[(a - 1) * k + kj] = temp;
                    }
                }
            }
            gettimeofday(&TP2, NULL);
            if (i == n - 1) printf("Sort Using time : %f ms\n", (TP2.tv_sec - TP1.tv_sec) * 1000.0 + (TP2.tv_usec - TP1.tv_usec) / 1000.0);
        }
        return;
    }

    // Recursively call Combination() to combine pivots
    int i;
    for (i = pivots[ki - 1] + 1; i < n; i++) {
        pivots[ki] = i;// p1 start from　 0　pivots里面所含的是索引
        //struct timeval TP1, TP2;
        //gettimeofday(&TP1, NULL);
        Combination(ki + 1, k, n, dim, M, coord, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
        //gettimeofday(&TP2, NULL);

        /** Iteration Log : pivots computed, best pivots, max distance sum, min distance sum pivots, min distance sum
        *** You can delete the logging code. **/
        if (ki == k - 2) {
            printf("Using time : %f ms\n", (TP2.tv_sec - TP1.tv_sec) * 1000.0 + (TP2.tv_usec - TP1.tv_usec) / 1000.0);
            int kj;
            printf("pivots: ");
            for (kj = 0; kj < k; kj++) {
                printf("%d \t", pivots[kj]);
            }
            printf("MaxPivots: ");
            for (kj = 0; kj < k; kj++) {
                printf("%d \t", maxDisSumPivots[kj]);
            }
            printf("maxDistanceSum: %lf\t", maxDistanceSum[0]);
            printf("MaxPivots: ");
            for (kj = 0; kj < k; kj++) {
                printf("%d \t", minDisSumPivots[kj]);
            }
            printf("minDistanceSum:%lf\n", minDistanceSum[0]);
        }
    }
}

int main(int argc, char* argv[]) {
    // filename : input file namespace
    char* filename = (char*)"uniformvector-2dim-5h.txt";
    if (argc == 2) {
        filename = argv[1];
    }
    else if (argc != 1) {
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
    if (file == NULL) {
        printf("%s file not found.\n", filename);
        return -1;
    }
    fscanf(file, "%d", &dim);//2
    fscanf(file, "%d", &n);//500
    fscanf(file, "%d", &k);//2
    printf("dim = %d, n = %d, k = %d\n", dim, n, k);

    // Read Data
    double* coord = (double*)malloc(sizeof(double) * dim * n);//500*2
    int i;
    for (i = 0; i < n; i++) {
        int j;
        for (j = 0; j < dim; j++) {
            fscanf(file, "%lf", &coord[i * dim + j]);//读取所有矩阵数据
        }
    }
    fclose(file);

    // Start timing
    struct timeval start;
    gettimeofday(&start, NULL);

    // maxDistanceSum : the largest M distance sum
    // minDistanceSum : the smallest M distance sum
    double* maxDistanceSum = (double*)malloc(sizeof(double) * (M + 1));//1001的空间,考虑存储空间,double 64
    double* minDistanceSum = (double*)malloc(sizeof(double) * (M + 1));
    for (i = 0; i < M; i++) {
        maxDistanceSum[i] = 0;
        minDistanceSum[i] = __DBL_MAX__;
    }

    // maxDisSumPivots : the top M pivots combinations
    // minDisSumPivots : the bottom M pivots combinations
    int* maxDisSumPivots = (int*)malloc(sizeof(int) * k * (M + 1));//2*1001
    int* minDisSumPivots = (int*)malloc(sizeof(int) * k * (M + 1));
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k; ki++) {
            maxDisSumPivots[i * k + ki] = 0;//创建了2*1001的为0的矩阵
            minDisSumPivots[i * k + ki] = 0;
        }
    }

    // temp : indexes of pivots with dummy array head
    int* temp = (int*)malloc(sizeof(int) * (k + 1));//pivots矩阵合集
    temp[0] = -1;

    //P2PDist : 存储点集里任意两点的距离的矩阵
    double* P2PDist = (double*)malloc(sizeof(double) * (n * n));
    SumP2PDistance(P2PDist, coord, n, dim);
    // Main loop. Combine different pivots with recursive function and evaluate them. Complexity : O( n^(k+2) )
    Combination(0, k, n, dim, M, coord, &temp[1], maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
    free(P2PDist);
    // End timing
    struct timeval end;
    gettimeofday(&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0);

    // Store the result
    FILE* out = fopen("result.txt", "w");
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", maxDisSumPivots[i * k + ki]);
        }
        fprintf(out, "%d\n", maxDisSumPivots[i * k + k - 1]);
    }
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", minDisSumPivots[i * k + ki]);
        }
        fprintf(out, "%d\n", minDisSumPivots[i * k + k - 1]);
    }
    fclose(out);

    // Log
    int ki;
    printf("max : ");
    for (ki = 0; ki < k; ki++) {
        printf("%d ", maxDisSumPivots[ki]);
    }
    printf("%lf\n", maxDistanceSum[0]);
    printf("min : ");
    for (ki = 0; ki < k; ki++) {
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
