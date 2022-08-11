#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<omp.h>
#include<unistd.h>

int main(){
    int i,j;
    omp_set_num_threads(8);
    omp_set_nested(8);
    for(i=0;i<8;i++){
        for(j=0;j<2;j++){
            printf("i=%d, j=%d, num_of_threads=%d\n",i,j,omp_get_thread_num);
        }
    }

}