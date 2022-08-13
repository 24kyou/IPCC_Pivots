#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<sys/time.h>
const int core_num = 32;
const int n = 500;
const int k = 2;

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
    printf("%d\n",loop);
    printf("cnt: %d, pivots:",cnt);
    for (j = 0; j < k; j++){
        start_pivots[cnt * k + j] = j;
        printf("%d ",start_pivots[cnt * k + j]);
    }
    printf("\n");
    OuterLoops++; cnt++;
    /*the loop is to caculate every start pivots[] which tend to start the main parallel*/
    for (i = k - 1; i >= 0; i--) {
        if (start_pivots[i] < i + n - k) {	//the last pivots[] num < n-1
            start_pivots[i]++;
            for (ki = i + 1; ki < k; ki++) start_pivots[ki] = start_pivots[ki - 1] + 1;
            i = k;
            if (OuterLoops == loop) {
                printf("cnt: %d, pivots:",cnt);
                for (j = 0; j < k; j++) {
                    start_pivots[cnt * k + j] = start_pivots[j];
                    printf("%d ", start_pivots[cnt * k + j]);
                }
                printf("\n");
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

int main(){
    int* start_pivots = (int*)malloc(sizeof(int)*k*core_num);
    SumStart_pivots(n,k,start_pivots);
}