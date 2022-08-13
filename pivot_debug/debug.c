#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<omp.h>

#define M 100000
#define N 25000

int main(int argc, char* argv[]){
    int i;
    printf("%d %d\n",M,N);
    int* B = (int*)malloc(sizeof(int)*M);
    //int* A = (int*)malloc(sizeof(int)*N);
    for(i=0;i<M;i++) B[i]=i;
    //for(i=0;i<N;i++) A[i]=0;
    i=0;
    for(i=0;i<M;i++) printf("%d ",B[i]);
    printf("\n");
    #pragma omp parallel num_threads(4)
    {
        int tid = omp_get_thread_num();
        //printf("0x%0x\n",&A);
        int *A = B+(tid*N);
        for(i=0;i<N;i++){
            A[i] = (4-tid)*N-i;
        }

        // for(i=0;i<N;i++){
        //     B[tid*N+i]=A[i];
        // }
        //改变printf的循环顺序会打乱B[]的读写操作
        // #pragma omp barrier
        // #pragma omp critical
        // {
        // for(i=0;i<N;i++){
        //     printf("%d %d\n",tid,A[i]);
        // }
        // }
        
    }
    printf("\n");
    for(i=0;i<M;i++){
        printf("%d ",B[i]);
        if((B[i]-1)%N==0) printf("\n");
    }
    printf("\n");
    free(B);
    return 0;
}																																														
