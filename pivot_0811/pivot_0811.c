#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<omp.h>

/*This fuction is to caculate every 2 point's distance and store it in the matrix P2PDist*/
void SumP2PDistance(double* P2PDist, double* coord, int n, int dim) {
    int i,ki;
#pragma omp parallel for private(i,ki) num_threads(32) collapse(2)
    for (i = 0; i < n; i++){
        for (ki = 0; ki < n; ki++) {
            double distance = 0;
            int j;
            for (j = 0; j < dim; j++) {
                distance += pow(coord[ki * dim + j] - coord[i * dim + j], 2);
            }
            //P2PDist[] Content: Distance[i,ki]
            P2PDist[i * n + ki] = sqrt(distance);
        }
    }
    return;
}

static unsigned long long OuterLoops = 0;
static unsigned long long one_loop = 1;
void SumStart_pivots(const int n, const int k, int* start_pivots) {
    int cnt = 0;
    unsigned long long loop;
    int i, j, ki;
    //struct timeval TP1, TP2;
    //gettimeofday(&TP1, NULL);
    for (i = 1; i <= k; i++) {
        one_loop *= (n - i + 1);
        one_loop /= i;
    }
    one_loop /= 32; //default core|threads number: 32;
    loop = one_loop;
    //printf("cnt: %d, pivots:",cnt);
    for (j = 0; j < k; j++) start_pivots[cnt * k + j] = j;
    //printf("\n");
    OuterLoops++; cnt++;
    /*the loop is to caculate every start pivots[] which tend to start the main parallel*/
    for (i = k - 1; i >= 0; i--) {
        if (start_pivots[i] < i + n - k) {	//the last pivots[] num < n-1
            start_pivots[i]++;
            for (ki = i + 1; ki < k; ki++) start_pivots[ki] = start_pivots[ki - 1] + 1;
            i = k;
            if (OuterLoops == loop) {
                //printf("cnt: %d, pivots:",cnt);
                for (j = 0; j < k; j++) {
                    //printf("%d ", pivots[j]);
                    start_pivots[cnt * k + j] = start_pivots[j];
                }
                //printf("\n");
                cnt++;
                if (cnt == 32) break;
                loop += one_loop;
            }
            OuterLoops++;
        }
    }
    for (j = 0; j < k; j++){
        start_pivots[j] = j;
        start_pivots[cnt * k + j] = n-j-k;
    }
    //gettimeofday(&TP2, NULL);
    //printf("Using time: %f ms\n", (TP2.tv_sec - TP1.tv_sec) * 1000.0 + (TP2.tv_usec - TP1.tv_usec) / 1000.0);
    return;
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


void caculate_sort(unsigned long long* inner_cnt, const int k, const int n, const int M, const int* pivots,
    double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots, const double* P2PDist){
    int i;
    double chebyshevSum = 0;
//#pragma omp parallel for private(i) reduction(+:chebyshevSum) num_threads(32) schedule(dynamic)
    for (i=0; i<n; i++){//using collapse  is slower
        int j;
	    for (j = i+1; j<n; j++) {   //openmp reduction sum
            double chebyshev = 0;
            int ki;
            for (ki = 0; ki < k; ki++) {    //  MAX(  |d(xi,pki)-d(xj,pki)| )   openmp reduction max
                double dis = fabs(P2PDist[i * n + pivots[ki]] - P2PDist[j * n + pivots[ki]]);
                chebyshev = dis > chebyshev ? dis : chebyshev;
            }
            chebyshevSum += chebyshev;
        }
    }
    chebyshevSum *= 2;
    //heap sort Time Complexity : O(M+(C_n^k-M)log^M)
    if(inner_cnt[0]<M-1){//Build heap, Time Complexity : O(M)
        maxDistanceSum[inner_cnt[0]] = chebyshevSum;
        minDistanceSum[inner_cnt[0]] = chebyshevSum;
        int kj;
        for (kj = 0; kj < k; kj++) {
            maxDisSumPivots[inner_cnt[0] * k + kj] = pivots[kj];
            minDisSumPivots[inner_cnt[0] * k + kj] = pivots[kj];
        }
    }else if(inner_cnt[0]==M-1){ //Adjust heap, Time Complexity : O(log^M)
        maxDistanceSum[inner_cnt[0]] = chebyshevSum;
        minDistanceSum[inner_cnt[0]] = chebyshevSum;
        int kj;
        for (kj = 0; kj < k; kj++) {
            maxDisSumPivots[inner_cnt[0] * k + kj] = pivots[kj];
            minDisSumPivots[inner_cnt[0] * k + kj] = pivots[kj];
        }
        int ki;
        for (ki = M / 2 - 1; ki >= 0; ki--) {
            make_minheap(maxDistanceSum, maxDisSumPivots, M-1, ki,k);
            make_maxheap(minDistanceSum, minDisSumPivots, M-1, ki, k);
        }
    }else{
        //Adjust heap, Time Complexity : O(log^M)
        int kj;
        if (chebyshevSum > maxDistanceSum[0]) {
            maxDistanceSum[0] = chebyshevSum;
            for (kj = 0; kj < k; kj++) maxDisSumPivots[kj] = pivots[kj];
            make_minheap(maxDistanceSum, maxDisSumPivots, M-1, 0, k);
        }
        if (chebyshevSum < minDistanceSum[0]) {
            minDistanceSum[0] = chebyshevSum;
            for (kj = 0; kj < k; kj++) minDisSumPivots[kj] = pivots[kj];
            make_maxheap(minDistanceSum, minDisSumPivots, M-1, 0, k);
        } 
    }
    inner_cnt[0]++;
}


void Combination(unsigned long long* inner_cnt, const int k, const int n, const int M, int* pivots, const int tid,
    double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots, const double* P2PDist) {
    int i, j;
    caculate_sort(inner_cnt, k, n, M, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
    for (i = k - 1; i >= 0; i--) {
        if (pivots[i] < i + n - k) {	//the least last maxDistanceSumay num < n-1
            pivots[i]++;
            for (j=i+1; j<k; j++) pivots[j] = pivots[j - 1] + 1;
            i = k;
            caculate_sort(inner_cnt, k, n, M, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
            if(tid==31) continue;   // last threads need to caculate mod part
            if(inner_cnt[0]==one_loop) break;
        }
    }
    return;
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
    int i,j;
    for (i = 0; i < n; i++) {
        int j;
        for (j = 0; j < dim; j++) {
            fscanf(file, "%lf", &coord[i * dim + j]);//��ȡ���о�������
        }
    }
    fclose(file);
    // Start timing
    struct timeval start;
    gettimeofday(&start, NULL);
    //these are shared param
    double* All_maxDistanceSum = (double*)malloc(sizeof(double) * M * 32);
    double* All_minDistanceSum = (double*)malloc(sizeof(double) * M * 32);
    int* All_maxDisSumPivots = (int*)malloc(sizeof(int) * k * M * 32);
    int* All_minDisSumPivots = (int*)malloc(sizeof(int) * k * M * 32);
    for (i = 0; i < M*32; i++) {
            All_maxDistanceSum[i] = 0;
            All_minDistanceSum[i] = __DBL_MAX__;
            for (j = 0; j < k; j++) {
                All_maxDisSumPivots[i * k + j] = 0;
                All_minDisSumPivots[i * k + j] = 0;
            }
    }
    //P2PDist : store every two point's distance
    double* P2PDist = (double*)malloc(sizeof(double) * (n * n));
    SumP2PDistance(P2PDist, coord, n, dim);
    //StartPivots: store the specific indexs of pivots to start with these threads  default:32core
    int* StartPivots = (int*)malloc(sizeof(int) * 32 * k);
    SumStart_pivots(n, k, StartPivots);
    // Main loop. Combine different pivots with recursive function and evaluate them. Complexity : O( n^(k+2) )  if(one_loop>...)
    //parallel
    #pragma omp parallel num_threads(32) 
    {
        int tid = omp_get_thread_num();
        int* pivots = StartPivots+tid*k;
        //StartPivots is a shared param    //broadcast the start pivots to every threads
        double* maxDistanceSum = All_maxDistanceSum+tid*M;
        double* minDistanceSum = All_minDistanceSum+tid*M;
        int* maxDisSumPivots = All_maxDisSumPivots+tid*M*k;
        int* minDisSumPivots = All_minDisSumPivots+tid*M*k;
        unsigned long long * inner_cnt =  (unsigned long long*)malloc(sizeof(unsigned long long));
        inner_cnt[0] = 0;

        caculate_sort(inner_cnt, k, n, M, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
        for (i = k - 1; i >= 0; i--) {
            if (pivots[i] < i + n - k) {	//the least last maxDistanceSumay num < n-1
                pivots[i]++;
                for (j=i+1; j<k; j++) pivots[j] = pivots[j - 1] + 1;
                i = k;
                #pragma omp critical
                {
                    caculate_sort(inner_cnt, k, n, M, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots, P2PDist);
                }
                
                if(tid==31) continue;   // last threads need to caculate mod part
                if(inner_cnt[0]==one_loop) break;
            }
        }
    }
    // // sort
    // int a;
    // for(i=0;i<32*M;i++){
    //     for (a = M*32-1; a > i; a--){
    //         if (All_maxDistanceSum[a] > All_maxDistanceSum[a - 1]){
    //             double temp = All_maxDistanceSum[a];
    //             All_maxDistanceSum[a] = All_maxDistanceSum[a - 1];
    //             All_maxDistanceSum[a - 1] = temp;
    //             int kj;
    //             for (kj = 0; kj < k; kj++){
    //                 int temp = All_maxDisSumPivots[a * k + kj];
    //                 All_maxDisSumPivots[a * k + kj] = All_maxDisSumPivots[(a - 1) * k + kj];
    //                 All_maxDisSumPivots[(a - 1) * k + kj] = temp;
    //             }
    //         }
    //         if (All_minDistanceSum[a] < All_minDistanceSum[a - 1]){
    //             double temp = All_minDistanceSum[a];
    //             All_minDistanceSum[a] = All_minDistanceSum[a - 1];
    //             All_minDistanceSum[a - 1] = temp;
    //             int kj;
    //             for (kj = 0; kj < k; kj++){
    //                 int temp = All_minDistanceSum[a * k + kj];
    //                 All_minDistanceSum[a * k + kj] = All_minDistanceSum[(a - 1) * k + kj];
    //                 All_minDistanceSum[(a - 1) * k + kj] = temp;
    //             }
    //         }
    //     }
    // }
    
    double* maxDistanceSum = All_maxDistanceSum;
    double* minDistanceSum = All_minDistanceSum;
    int* maxDisSumPivots = All_maxDisSumPivots;
    int* minDisSumPivots = All_maxDisSumPivots;
    // //heap sort again
    // //1st: build heap 
    // for (i = M / 2 - 1; i >= 0; i--) {
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
    // double* maxDistanceSum = All_maxDistanceSum;
    // double* minDistanceSum = All_minDistanceSum;
    // int* maxDisSumPivots = All_maxDisSumPivots;
    // int* minDisSumPivots = All_maxDisSumPivots;
    // End timing
    struct timeval end;
    gettimeofday(&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0);

    //free
    //free(P2PDist); free(All_maxDistanceSum); free(All_minDistanceSum); free(All_maxDisSumPivots); free(All_minDisSumPivots);

    // Store the result
    FILE* out = fopen("result.txt", "w");
    for (i = 0; i < M*32; i++) {
        int ki;
        for (ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", maxDisSumPivots[i * k + ki]);
        }
        fprintf(out, "%d ", maxDisSumPivots[i * k + k - 1]);
        fprintf(out, "%lf\n", maxDistanceSum[i]);
    }
    for (i = 0; i < M*32; i++) {
        int ki;
        for (ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", minDisSumPivots[i * k + ki]);
        }
        fprintf(out, "%d ", minDisSumPivots[i * k + k - 1]);
        fprintf(out, "%lf\n", minDistanceSum[i]);
    }
    fclose(out);
    //free(maxDisSumPivots);free(minDisSumPivots);free(maxDistanceSum);free(minDistanceSum);free(coord);
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
