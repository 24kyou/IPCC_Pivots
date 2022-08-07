#include<stdio.h>
#include<omp.h>
#include<time.h>
#include<stdlib.h>

int Partition(int* data, int start, int end)   //划分数据
{
    int temp = data[start];   //以第一个元素为基准
    while (start < end) {
        while (start < end && data[end] >= temp)end--;   //找到第一个比基准小的数
        data[start] = data[end];
        while (start < end && data[start] <= temp)start++;    //找到第一个比基准大的数
        data[end] = data[start];
    }
    data[start] = temp;   //以基准作为分界线
    return start;
}

void quickSort(int* data, int start, int end)  //并行快排
{
    if (start < end) {
        int pos = Partition(data, start, end);
        #pragma omp parallel sections    //设置并行区域
        {
            #pragma omp section          //该区域对前部分数据进行排序
            quickSort(data, start, pos - 1);
            #pragma omp section          //该区域对后部分数据进行排序
            quickSort(data, pos + 1, end);
        }
    }
}

int main(int argc, char* argv[])
{
    int n = atoi(argv[2]), i;   //线程数
    int size = atoi(argv[1]);   //数据大小
    int* num = (int*)malloc(sizeof(int) * size);

    double starttime = omp_get_wtime();
    srand(time(NULL) + rand());   //生成随机数组
    for (i = 0; i < size; i++)
        num[i] = rand();
    omp_set_num_threads(n);   //设置线程数
    quickSort(num, 0, size - 1);   //并行快排
    double endtime = omp_get_wtime();

    for (i = 0; i < 10 && i<size; i++)//输出前十个元素
        printf("%d ", num[i]);
    printf("\n并行时间：%lfs\n", endtime - starttime);
    return 0;
}
