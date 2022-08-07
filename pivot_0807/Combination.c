#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<sys/time.h>
#define  MAX 100

typedef void OutProc(int[], int);

// output
void  OutputCombination(int  ary[], int  n) {
    static int count = 0;
    int i;
    printf("%05d : ", ++count);
    for (i = 0; i < n; i++) {
        printf("%d ", ary[i]);
    }
    printf(" \n");
}

void TestCombination(int n, int r, OutProc proc) {
    int ary[MAX];
    int i, k;
    struct timeval TP1,TP2;
    gettimeofday(&TP1,NULL);
    for (i = 0; i < r; i++) ary[i] = i;
    proc(ary, r);
    for (i = r - 1; i >= 0; i--) {
        if (ary[i] < i + n - r) {	//the least last array num < n-1
            ary[i]++;
            for (k = i + 1; k < r; k++) ary[k] = ary[k - 1] + 1;
            i = r;
            proc(ary, r);
        }
    }
    gettimeofday(&TP2,NULL);
    printf("Using time: %f ms\n",(TP2.tv_sec-TP1.tv_sec)*1000.0+(TP2.tv_usec-TP1.tv_usec)/1000.0);
}

// main
int  main() {
    TestCombination(20, 4, OutputCombination);
    return 0;
}
