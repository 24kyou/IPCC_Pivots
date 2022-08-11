#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<omp.h>

int main(int argc, char* argv[]){
    int i;
    int* B = (int*)malloc(sizeof(int)*8);
    for(i=0;i<8;i++) B[i]=i;
    for(i=0;i<8;i++) printf("%d ",B[i]);
    printf("\n");
    #pragma omp parallel num_threads(4) 
    {
        int* A = (int*)malloc(sizeof(int)*2);
        int tid = omp_get_thread_num();
        for(i=0;i<2;i++){
            A[i] = B[tid*2+i];
        }
        printf("threads: %d, &A: 0x%lx, Read A[0] A[1] as %d %d\n",tid,&A,A[0],A[1]);
        //B[tid*2] = A[0];
        //B[tid*2+1] = A[1];
        //printf("threads: %d, &B: 0x%lx, Read B[tid*2] B[tid*2+1] as %d %d\n",tid,&B,B[tid*2],B[tid*2+1]);
        free(A);
    }
    for(i=0;i<8;i++) printf("%d ",B[i]);
    printf("\n");
    printf("over\n");
    free(B);
    return 0;
}																																														
