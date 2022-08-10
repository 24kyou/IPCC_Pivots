/*******
    run.sh:
    #!/bin/bash
    #SBATCH --job-name=Pivot0809
    #SBATCH --partition=IPCC //如果在学校超算为 cu
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=32
    #SBATCH --error=%j.err
    #SBATCH --output=%j.out
    icpc pivot_0808.c -Ofast -lm -qopenmp -o pivot_0809
    ./pivot_0809

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
            //maxDistanceSumay Content: Distance[i,ki]
            P2PDist[i * n + ki] = sqrt(distance);
            //Testlog:
            //if(i==0) printf("point1:%d\tpoint2:%d\tdistance:%lf\n", i, ki, P2PDist[i * (n - 1) + ki - subnum]);
        }
    }
    return;
}

static unsigned long long count = 0;
void SumStart_pivots(const int n, const int k, int* start_pivots,int* pivots) {
    int cnt = 0;
    unsigned long long loop,one_loop=1;
    int i, j, ki;
    //struct timeval TP1, TP2;
    //gettimeofday(&TP1, NULL);

    for (i = 1; i <= k; i++) {
        one_loop *= (n - i + 1);
        one_loop /= i;
    }
    one_loop /= 32; //default core 32;
    /*这里很重要!!具体怎么分取决于你的数据输入量*/
    loop = one_loop;

    for (i = 0; i < k; i++) pivots[i] = i;
    for (j = 0; j < k; j++) {
        //printf("%d ", pivots[j]);
        start_pivots[cnt * k + j] = pivots[j];
    }
    //printf("\n");
    count++; cnt++;

    for (i = k - 1; i >= 0; i--) {
        if (pivots[i] < i + n - k) {	//the least last maxDistanceSumay num < n-1
            pivots[i]++;
            for (ki = i + 1; ki < k; ki++) pivots[ki] = pivots[ki - 1] + 1;
            i = n;
            if (count == loop) {
                //printf("%lld: %lld : ",cnt,count);
                //printf("%lld: ", cnt);
                for (j = 0; j < n; j++) {
                    //printf("%d ", ary[j]);
                    start_pivots[cnt * n + j] = pivots[j];
                }
                //printf("\n");
                cnt++;
                if (cnt == 64) break;
                loop += one_loop;
            }
            count++;
        }
    }
    return;
    //gettimeofday(&TP2, NULL);
    //printf("Using time: %f ms\n", (TP2.tv_sec - TP1.tv_sec) * 1000.0 + (TP2.tv_usec - TP1.tv_usec) / 1000.0);
}

void make_maxheap(double* minDistanceSum, int* minDisSumPivots, int end, int start,int k) {
    int fa = start;
    int child = fa * 2 + 1;
    while (child <= end) {
        if (child + 1 <= end && minDistanceSum[child] < minDistanceSum[child + 1]) child++;
        if (minDistanceSum[fa] > minDistanceSum[child]) return;
        else {
            double temp = minDistanceSum[child];
            minDistanceSum[child] = minDistanceSum[fa];
            minDistanceSum[fa] = temp;
            int kj;
            for (kj = 0; kj < k; kj++) {
                int temp = minDisSumPivots[child * k + kj];
                minDisSumPivots[child * k + kj] = minDisSumPivots[fa * k + kj];
                minDisSumPivots[fa * k + kj] = temp;
            }
            fa = child;
            child = fa * 2 + 1;
        }
    }
}

void make_minheap(double* maxDistanceSum, int* maxDisSumPivots, int end, int start, int k) {
    int fa = start;
    int child = fa * 2 + 1;
    while (child <= end) {
        if (child + 1 <= end && maxDistanceSum[child] > maxDistanceSum[child + 1]) child++;
        if (maxDistanceSum[fa] < maxDistanceSum[child]) return;
        else{
            double temp = maxDistanceSum[child];
            maxDistanceSum[child] = maxDistanceSum[fa];
            maxDistanceSum[fa] = temp;
            int kj;
            for (kj = 0; kj < k; kj++) {
                int temp = maxDisSumPivots[child * k + kj];
                maxDisSumPivots[child * k + kj] = maxDisSumPivots[fa * k + kj];
                maxDisSumPivots[fa * k + kj] = temp;
            }
            fa = child;
            child = fa * 2 + 1;
        }
    }
}

int cnt = 0;
int flag = 1;
void caculate_sort(const int k, const int n, const int M, int* pivots,
    double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots, double* P2PDist){
    int i;
    //caculate
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
    if (pivots[0] == 46 && pivots[1] == 77 && pivots[2] == 80 && pivots[3] == 92 && pivots[4] == 99) printf(" 46 77 80 92 99 : %lf\n", chebyshevSum);
    if (pivots[0] == 16 && pivots[1] == 35 && pivots[2] == 77 && pivots[3] == 80 && pivots[4] == 99) printf("16 35 77 80 99 : %lf\n", chebyshevSum);
    if (pivots[0] == 6 && pivots[1] == 31 && pivots[2] == 46 && pivots[3] == 77 && pivots[4] == 80) printf(" 6 31 46 77 80 : %lf\n", chebyshevSum);
 /*   46 77 80 92 99
    16 35 77 80 99
    6 31 46 77 80*/
    //sort
    if(cnt<M){
        maxDistanceSum[cnt] = chebyshevSum;
        minDistanceSum[cnt] = chebyshevSum;
        int kj;
        for (kj = 0; kj < k; kj++) {
            maxDisSumPivots[cnt * k + kj] = pivots[kj];
            minDisSumPivots[cnt * k + kj] = pivots[kj];
        }
        cnt++;
    }
    if(cnt==M && flag){
        int ki;
        //1st: build heap 
        for (ki = M / 2 - 1; ki >= 0; ki--) {
            make_minheap(maxDistanceSum, maxDisSumPivots, M-1, ki,k);
            make_maxheap(minDistanceSum, minDisSumPivots, M-1, ki, k);
        }
        flag--;
    }
    else{
        int kj;
#pragma omp sections
        {
#pragma omp section
            if (chebyshevSum > maxDistanceSum[0]) {
                maxDistanceSum[0] = chebyshevSum;
                for (kj = 0; kj < k; kj++) maxDisSumPivots[kj] = pivots[kj];
                make_minheap(maxDistanceSum, maxDisSumPivots, M, 0, k);
            }
#pragma omp section
            if (chebyshevSum < minDistanceSum[0]) {
                minDistanceSum[0] = chebyshevSum;
                for (kj = 0; kj < k; kj++) minDisSumPivots[kj] = pivots[kj];
                make_maxheap(minDistanceSum, minDisSumPivots, M, 0, k);

            }
        }
    }
}


void Combination(const int k, const int n, const int M, int* pivots,
    double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots, double* P2PDist) {
    int i, j;
    //第一次
    caculate_sort(k, n, M, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
    for (i = k - 1; i >= 0; i--) {
        if (pivots[i] < i + n - k) {	//the least last maxDistanceSumay num < n-1
            pivots[i]++;
            for (j=i+1; j<k; j++) pivots[j] = pivots[j - 1] + 1;
            i = k;
            
            caculate_sort(k, n, M, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
        }
    }
    //sort minheap
    for (i = M - 1; i >= 0; i--) {
        //swap(0,i)
        double temp = maxDistanceSum[i];
        maxDistanceSum[i] = maxDistanceSum[0];
        maxDistanceSum[0] = temp;
        int j;
        for (j = 0; j < k; j++) {
            int temp = maxDisSumPivots[i * k + j];
            maxDisSumPivots[i * k + j] = maxDisSumPivots[j];
            maxDisSumPivots[j] = temp;
        }
        make_minheap(maxDistanceSum, maxDisSumPivots, i - 1, 0, k);
    }
    //sort maxheap 
    for (i = M - 1; i >= 0; i--) {
        //swap(0,i)
        double temp = minDistanceSum[i];
        minDistanceSum[i] = minDistanceSum[0];
        minDistanceSum[0] = temp;
        int j;
        for (j = 0; j < k; j++) {
            int temp = minDisSumPivots[i * k + j];
            minDisSumPivots[i * k + j] = minDisSumPivots[j];
            minDisSumPivots[j] = temp;
        }
        make_maxheap(minDistanceSum, minDisSumPivots, i - 1, 0, k);
    }
}

int main(int argc, char* argv[]){
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
    fscanf(file, "%d", &dim);
    fscanf(file, "%d", &n);
    fscanf(file, "%d", &k);
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
    int* maxDisSumPivots = (int*)malloc(sizeof(int) * k * (M + 1));
    int* minDisSumPivots = (int*)malloc(sizeof(int) * k * (M + 1));
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k; ki++) {
            maxDisSumPivots[i * k + ki] = 0;
            minDisSumPivots[i * k + ki] = 0;
        }
    }

    // pivots : indexes of pivots with dummy maxDistanceSumay head
    int* pivots = (int*)malloc(sizeof(int) * k);
    for (i = 0; i < k; i++) pivots[i] = i;

    //P2PDist : store every two point's distance
    double* P2PDist = (double*)malloc(sizeof(double) * (n * n));
    SumP2PDistance(P2PDist, coord, n, dim);

    //StartPivots: store the specific indexs of pivots to start with these threads  default:32core
    int* StartPivots = (int*)malloc(sizeof(int) * 31 * k);
    //SumStart_pivots(n, k, StartPivots,pivots);
    for (i = 0; i < k; i++) pivots[i] = i;

    // Main loop. Combine different pivots with recursive function and evaluate them. Complexity : O( n^(k+2) )
    Combination(k, n, M, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);

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

    return 0;
}
