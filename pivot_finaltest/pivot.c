#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<omp.h>

#define OuterThreads 32

/*This fuction is to caculate every 2 point's distance and store it in the matrix P2PDist*/
void SumP2PDistance(double* P2PDist, double* coord, int n, int dim) {
    int i,ki;
//#pragma omp parallel for private(i,ki) num_threads(32) collapse(2)
    for (i = 0; i < n; i++){
        for (ki = 0; ki < n; ki++) {
            double distance = 0;
            int j;
            for (j = 0; j < dim; j++) {
                distance += pow(coord[ki * dim + j] - coord[i * dim + j], 2);
            }
            //P2PDist[] Content: Distance[i,ki]
            P2PDist[i * n + ki] = sqrt(distance);
            //Testlog:
            //if(i==0) printf("point1:%d\tpoint2:%d\tdistance:%lf\n", i, ki, P2PDist[i * (n - 1) + ki - subnum]);
        }
    }
    return;
}
/*This function is to assure every OuterThreads's StartPivots[] and caculate*/
void SumStart_pivots(const int n, const int k, int* StartPivots, int*OneLoopSum, int*OuterLoopSum) {
    int cnt = 0;
    OuterLoopSum[0]=1; OneLoopSum[0]=0;
    int i, j, ki;
    //struct timeval TP1, TP2;
    //gettimeofday(&TP1, NULL);
    for (i = 1; i <= k; i++) {
        OuterLoopSum[0] *= (n - i + 1);
        OuterLoopSum[0] /= i;
    }
    OneLoopSum[0] = OuterLoopSum[0]/OuterThreads;
    //printf("%d %d",OuterLoopSum[0],OneLoopSum[0]);
    int InnerLoopSum=0;
    //printf("cnt: %d, pivots:",cnt);
    for (i = 0; i < k; i++) {
        StartPivots[i] = i;
    }
    InnerLoopSum++; cnt++;
    /*the loop is to caculate every start pivots[] which tend to start the main parallel*/
    for (i = k - 1; i >= 0; i--) {
        if (StartPivots[i] < i + n - k) {	//the last pivots[] num < n-1
            StartPivots[i]++;
            for (ki = i + 1; ki < k; ki++) StartPivots[ki] = StartPivots[ki - 1] + 1;
            i = k;
            if (InnerLoopSum == OneLoopSum[0]) {
                for (j = 0; j < k; j++) {
                    StartPivots[cnt * k + j] = StartPivots[j];
                }
                //printf("\n");
                cnt++;
                if (cnt == OuterThreads) break;
                OneLoopSum[0] += OuterLoopSum[0]/OuterThreads;
            }
            InnerLoopSum++;
        }
    }
    OneLoopSum[0] = OuterLoopSum[0]/OuterThreads;
    for (i = 0; i < k; i++) {
        StartPivots[i] = i;
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

void caculate_sort(int* innerCnt, const int k, const int n, const int M, int* pivots,
    double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots, double* P2PDist){
    int i;
    //caculate
    double chebyshevSum = 0;
//#pragma omp parallel for private(i) reduction(+:chebyshevSum) num_threads(32)
    for (i = 0; i < n; i++) {
        int j;
	for (j = i+1; j < n; j++) {   //openmp reduction sum
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
    //sort
    if(innerCnt[0]<M-1){
        maxDistanceSum[innerCnt[0]] = chebyshevSum;
        minDistanceSum[innerCnt[0]] = chebyshevSum;
        int kj;
        for (kj = 0; kj < k; kj++) {
            maxDisSumPivots[innerCnt[0] * k + kj] = pivots[kj];
            minDisSumPivots[innerCnt[0] * k + kj] = pivots[kj];
        }
    }
    if(innerCnt[0]==M-1){
        maxDistanceSum[innerCnt[0]] = chebyshevSum;
        minDistanceSum[innerCnt[0]] = chebyshevSum;
        int kj;
        for (kj = 0; kj < k; kj++) {
            maxDisSumPivots[innerCnt[0] * k + kj] = pivots[kj];
            minDisSumPivots[innerCnt[0] * k + kj] = pivots[kj];
        }
        int ki;
        //1st: build heap 
        for (ki = M / 2 - 1; ki >= 0; ki--) {
            make_minheap(maxDistanceSum, maxDisSumPivots, M-1, ki,k);
            make_maxheap(minDistanceSum, minDisSumPivots, M-1, ki, k);
        }
    }
    else{
        int kj;
//#pragma omp sections
        //{
//#pragma omp section
            if (chebyshevSum > maxDistanceSum[0]) {
                maxDistanceSum[0] = chebyshevSum;
                for (kj = 0; kj < k; kj++) maxDisSumPivots[kj] = pivots[kj];
                make_minheap(maxDistanceSum, maxDisSumPivots, M-1, 0, k);
            }
//#pragma omp section
            if (chebyshevSum < minDistanceSum[0]) {
                minDistanceSum[0] = chebyshevSum;
                for (kj = 0; kj < k; kj++) minDisSumPivots[kj] = pivots[kj];
                make_maxheap(minDistanceSum, minDisSumPivots, M-1, 0, k);

            }
        //}
    }
    innerCnt[0]++;
}


void Combination(int* innerCnt, const int k, const int n, const int M, int* pivots,
    double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots, double* P2PDist) {
    int i, j;
    innerCnt[0]=0;
    caculate_sort(innerCnt,k, n, M, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
    for (i = k - 1; i >= 0; i--) {
        if (pivots[i] < i + n - k) {	//the least last maxDistanceSumay num < n-1
            pivots[i]++;
            for (j=i+1; j<k; j++) pivots[j] = pivots[j - 1] + 1;
            i = k;
            caculate_sort(innerCnt,k, n, M, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
        }
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
    double* coord = (double*)malloc(sizeof(double) * dim * n);
    int i,j;
    for (i = 0; i < n; i++) {
        int j;
        for (j = 0; j < dim; j++) {
            fscanf(file, "%lf", &coord[i * dim + j]);
        }
    }
    fclose(file);

    // Start timing
    struct timeval start;
    gettimeofday(&start, NULL);

    // maxDistanceSum : the largest M distance sum
    // minDistanceSum : the smallest M distance sum
    // maxDisSumPivots : the top M pivots combinations
    // minDisSumPivots : the bottom M pivots combinations
    double* All_maxDistanceSum = (double*)malloc(sizeof(double) * M * OuterThreads);
    double* All_minDistanceSum = (double*)malloc(sizeof(double) * M * OuterThreads);
    int* All_maxDisSumPivots = (int*)malloc(sizeof(int) * k * M * OuterThreads);
    int* All_minDisSumPivots = (int*)malloc(sizeof(int) * k * M * OuterThreads);
    for (i = 0; i < M*OuterThreads; i++) {
        All_maxDistanceSum[i] = 0;
        All_minDistanceSum[i] = 0;
        int ki;
        for (ki = 0; ki < k; ki++) {
            All_maxDisSumPivots[i * k + ki] = 0;
            All_minDisSumPivots[i * k + ki] = 0;
        }
    }
    //P2PDist : store every two point's distance
    double* P2PDist = (double*)malloc(sizeof(double) * (n * n));
    SumP2PDistance(P2PDist, coord, n, dim);
    //StartPivots: store the specific indexs of pivots to start with these threads  default:32core
    int* StartPivots = (int*)malloc(sizeof(int) * 32 * k);    
    int* OneLoopSum = (int*)malloc(sizeof(int));
    int* OuterLoopSum = (int*)malloc(sizeof(int));
    SumStart_pivots(n, k, StartPivots,OneLoopSum,OuterLoopSum);
    // Main loop. Combine different pivots with recursive function and evaluate them. Complexity : O( n^(k+2) )
    #pragma omp parallel num_threads(OuterThreads) shared(StartPivots, All_maxDistanceSum, All_minDistanceSum, All_maxDisSumPivots, All_minDisSumPivots)
    {
        int tid = omp_get_thread_num();
        int* innerCnt = (int*)malloc(sizeof(int));
        innerCnt[0]=0;
        int* pivots = (int*)malloc(sizeof(int)*k);
        #pragma omp critical
        {
            for(i=0;i<k;i++) pivots[i] = StartPivots[tid*k+i];
        }
        double* maxDistanceSum = (double*)malloc(sizeof(double) * M);
        double* minDistanceSum = (double*)malloc(sizeof(double) * M);
        int* maxDisSumPivots = (int*)malloc(sizeof(int) * k * M);
        int* minDisSumPivots = (int*)malloc(sizeof(int) * k * M);
        for (i = 0; i < M; i++) {
            maxDistanceSum[i] = 0;
            minDistanceSum[i] = 0;
            int ki;
            for (ki = 0; ki < k; ki++) {
                maxDisSumPivots[i * k + ki] = 0;
                minDisSumPivots[i * k + ki] = 0;
            }
        }
        Combination(innerCnt, k, n, M, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
        #pragma omp critical
        {
        int a = tid*M;
        for(i=0;i<M;i++){
            All_maxDistanceSum[a+i] = maxDistanceSum[i];
            All_minDistanceSum[a+i] = minDistanceSum[i];
            int ki;
            for(ki=0;ki<k;ki++){
                All_maxDisSumPivots[(a+i)*k+ki] = maxDisSumPivots[i*k+ki];
                All_maxDisSumPivots[(a+i)*k+ki] = maxDisSumPivots[i*k+ki];
            }
        }
        }
    }
    // //heap sort
    // //1st: build heap 
    // for (i = M/2 - 1; i >= 0; i--) {
    //     make_minheap(All_maxDistanceSum, All_maxDisSumPivots, M-1, i, k);
    //     make_maxheap(All_minDistanceSum, All_minDisSumPivots, M-1, i, k);
    // }
    // //2nd: adjust heap
    // for(i=M;i<M*32;i++){
    //     if (All_maxDistanceSum[i] > All_maxDistanceSum[0]) {
    //         All_maxDistanceSum[0] = All_maxDistanceSum[i];
    //         for (j = 0; j < k; j++) All_maxDistanceSum[j] = All_maxDistanceSum[i*k+j];
    //         make_minheap(All_maxDistanceSum, All_maxDisSumPivots, M-1, 0, k);
    //     }
    //     if (All_minDistanceSum[i] < All_minDistanceSum[0]) {
    //         All_minDistanceSum[0] = All_minDistanceSum[i];
    //         for (j = 0; j < k; j++) All_minDistanceSum[j] = All_minDistanceSum[i*k+j];
    //         make_maxheap(All_minDistanceSum, All_minDisSumPivots, M-1, 0, k);
    //     }
    // }
    // //sort
    // for (i = M - 1; i >= 0; i--) {
    //     //swap(0,i)
    //     double temp = All_maxDistanceSum[i];
    //     All_maxDistanceSum[i] = All_maxDistanceSum[0];
    //     All_maxDistanceSum[0] = temp;
    //     temp = All_minDistanceSum[i];
    //     All_minDistanceSum[i] = All_minDistanceSum[0];
    //     All_minDistanceSum[0] = temp;
    //     int j;
    //     for (j = 0; j < k; j++) {
    //         int temp = All_maxDisSumPivots[i * k + j];
    //         All_maxDisSumPivots[i * k + j] = All_maxDisSumPivots[j];
    //         All_maxDisSumPivots[j] = temp;
    //         temp = All_minDisSumPivots[i * k + j];
    //         All_minDisSumPivots[i * k + j] = All_minDisSumPivots[j];
    //         All_minDisSumPivots[j] = temp;
    //     }
    //     make_minheap(All_maxDistanceSum, All_maxDisSumPivots, i - 1, 0, k);
    //     make_maxheap(All_minDistanceSum, All_minDisSumPivots, i - 1, 0, k);
    // }
    free(P2PDist);
    // End timing
    struct timeval end;
    gettimeofday(&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0);
    // Store the result
    FILE* out = fopen("result.txt", "w");
    for (i = 0; i < M*32; i++) {
        int ki;
        for (ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", All_maxDisSumPivots[i * k + ki]);
        }
        fprintf(out, "%d", All_maxDisSumPivots[i * k + k - 1]);
        fprintf(out, "%lf\n", All_maxDistanceSum[i]);
    }
    for (i = 0; i < M*32; i++) {
        int ki;
        for (ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", All_minDisSumPivots[i * k + ki]);
        }
        fprintf(out, "%d", All_minDisSumPivots[i * k + k - 1]);
        fprintf(out, "%d\n", All_minDistanceSum[i]);
    }
    fclose(out);
    // Log
    int ki;
    printf("max : ");
    for (ki = 0; ki < k; ki++) {
        printf("%d ", All_maxDisSumPivots[ki]);
    }
    printf("%lf\n", All_maxDistanceSum[0]);
    printf("min : ");
    for (ki = 0; ki < k; ki++) {
        printf("%d ", All_minDisSumPivots[ki]);
    }
    printf("%lf\n", All_minDistanceSum[0]);
    return 0;
}																																														