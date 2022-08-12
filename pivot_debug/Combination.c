#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<sys/time.h>

static unsigned long long count = 0;

unsigned long long Combination(int n, int m){
    //C^n_m
    unsigned long long res=1;
    int i;    
    for(i=1;i<=n;i++){
        res*=(m-i+1);
        res/=i;
    }
    return res;
}

void TestCombination(int m, int n, int*start_pivots) {
    int cnt=0;
    unsigned long long one_loop= Combination(n, m)/64;
    int i, j, k;
    struct timeval TP1,TP2;
    gettimeofday(&TP1,NULL);

    for (i = 0; i < n; i++) ary[i] = i;
    printf("%lld: ",cnt);
    for(j=0;j<n;j++){
        printf("%d ",ary[j]);
	start_pivots[cnt*n+j] = ary[j];
    }
printf("\n");
    count++;
    cnt++;

    for (i = n - 1; i >= 0; i--) {
        if (ary[i] < i + m - n) {	//the least last array num < n-1
            ary[i]++;
            for (k = i + 1; k < n; k++) ary[k] = ary[k - 1] + 1;
            i = n;

            if(count==loop){
		//printf("%lld: %lld : ",cnt,count);
		printf("%lld: ",cnt);
		for(j=0;j<n;j++){
		    printf("%d ",ary[j]);
	            start_pivots[cnt*n+j] = ary[j];
		}
		printf("\n");
		cnt++;
		if(cnt==64) break;
		loop+=one_loop;
            }
            count++;
        }
    }
    gettimeofday(&TP2,NULL);
    printf("Using time: %f ms\n",(TP2.tv_sec-TP1.tv_sec)*1000.0+(TP2.tv_usec-TP1.tv_usec)/1000.0);
}

// main
int  main() {
    int n,m;//n<m
    printf("Input your n and m\n");
    scanf("%d %d",&n,&m);
    int* start_pivots = (int*)malloc(sizeof(int)*63*n);
    //unsigned long long res = Combination(n,m);
    //unsigned long long one_loop = res/64;
    //printf("All_loopsum:%lld AVE_loopsum:%lld\n",res,one_loop);
    TestCombination(m,n,start_pivots);
    int i,j; 
    //outputlog
    for(i=0;i<64;i++){
	    for(j=0;j<n;j++){
		    printf("%d ",start_pivots[i*n+j]);
	    }
	    printf("\n");
    }
    return 0;
}
