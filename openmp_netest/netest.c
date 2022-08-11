#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<omp.h>
#include<unistd.h>
#define  ull unsigned long long

ull Combination(int n,int k){
    int i;
    ull loop=1;
    for (i = 1; i <= k; i++) {
        loop *= (n - i + 1);
        loop /= i;
    }
    return loop;
}

int main(){
    //printf("%d",omp_get_thread_num);
    ull outer, inner;
    int outer_n,outer_k,inner_n,inner_k;

    printf("Input the outer Combination n and k\n");
    scanf("%d %d",&outer_n,&outer_k);
    printf("Input the inner Combination n and k\n");
    scanf("%d %d",&inner_n,&inner_k);

    outer =  Combination(outer_n,outer_k);
    inner = Combination(inner_n,inner_k);
    printf("%lld %lld\n",outer,inner);
    ull count =0;
    ull i,j;
    struct timeval tp1,tp2;
    struct timeval tp3,tp4;
    double timecnt1=0,timecnt2=0;
    
    gettimeofday(&tp3,NULL);
    omp_set_nested(32);
    #pragma omp parallel for private(j) reduction(+:count,timecnt1) num_threads(16) 
    for(i=0;i<outer;i++){
        gettimeofday(&tp1,NULL);
        #pragma omp parallel for reduction(+:count) num_threads(2) 
        for(j=0;j<inner;j++){
            count++;
            //sleep(10);
        }
        gettimeofday(&tp2,NULL);
        timecnt1+=(tp2.tv_sec - tp1.tv_sec) * 1000.0 + (tp2.tv_usec - tp1.tv_usec) / 1000.0;
        //printf("timecn1: %lfms\n",timecnt1);
        //printf("%d\n",count);
    }
    gettimeofday(&tp4,NULL);
    printf("timecn1: %lfms  timecnt2: %lfms\n",timecnt1,(tp4.tv_sec - tp3.tv_sec) * 1000.0 + (tp4.tv_usec - tp3.tv_usec) / 1000.0);
    printf("%lld\n",count);
}