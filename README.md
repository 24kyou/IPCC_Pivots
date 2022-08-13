# 2022_IPCC 初赛赛题优化报告总结



### 0707

学习Git的基本使用-->>版本迭代

在自己的Github上构建自己的repo

安装Visual Studio	增加mpi函数库环境

开始学习简单的利用mpi进行多线程编程

部分资料：

[(Windows系统下Visual studio 2022MPI 环境配置](https://blog.csdn.net/weixin_45965810/article/details/124900138)

BiliBili ：MPI系列

[【MPI系列】第一节-MPI并行编程技术-基本概念_哔哩哔哩_bilibili](https://www.bilibili.com/video/BV1Mg411P7Bt?spm_id_from=333.999.0.0)

[MPI 常用函数概述](https://blog.csdn.net/weixin_40255793/article/details/84201243)

**接下来的：**

MPI编程	并行程序与算法	OpenMP

### 0711

mpi库函数学习

### 0712

Makefile学习

### 0715

制定优化方案，

1.通过编译选项，加速库

2.对源程序的算法进行优化调整

### 0719

写完第一个mpidemo

### 0720

初赛题发布

尝试使用学校超算资源跑一遍

需要开始学习一些操作指令  ssh  slurm sbatch 也即超算的基本使用

使用scp在linux服务器之间进行文件传输

以本人的scp 从 本地上传到超算平台 举例：

```powershell
# ls
pivot.tar
# pwd
/root/Desktop
# scp /root/Desktop pivot.tar hr@10.160.195.152:/public/home/hr
/*输入密码即可*/
/*scp <文件目录> <文件名> <账户>@<Linux服务器的IP地址>：<此服务器中你要上传文件的所在目录>*/
```

简单看了一下试题，有初步的几个优化点；

1. 最开始的编译器选项优化

2. 循环合并 赋值循环计算循环都有比较大的合并点
3. 拆解递归 把那个函数的递归拆成普通循环
4. 对循环做openmp的并行化
5. 在这些基础上看看能不能对具体计算做一些优化

## 实验数据记录&优化过程记录

初始数据:

跑一次Combination大循环大概是1413ms 未优化

跑一次distanceSum大循环大概是2-4ms 未优化

**得到结果：**

**max:143 351 117761.647418**

**min:83 226 43769.849602**

## 赛题源代码（不经任何优化）耗时：440104.109ms

## 针对解题报告里的优化点3进行优化 得到耗时：183401.57ms

## 开O3优化 得到耗时：33002.64200ms

## 开O2优化 得到耗时：28000ms

### 0729:

到周一0801前做的三个点（希望不是给自己画饼...)

1. 重构Sum_distance函数 新建距离矩阵和计算
2. 将原先的堆排序改为并行优化的快速排序
3. 对大循环进行openmp优化（不一定）
4. 针对递归 思考是否有循环化的可能性（时间复杂度会很高）。

做完前两条就算胜利

如何快速检验两个文本文件是否相等：

```powershell
md5 file1
<md5_code> file1
md5 file2
<md5_code> file2
/*此后只需观察两个code是否相等即可*/
```

## 现在对解题思路的大致思考：

获得原题->进行题目解读->算法优化->并行优化->细节数据结构优化->结束

### 0730

针对点1的思考

第一次尝试：优化rebuildcoord矩阵 

第二次尝试：建立P2PDistance矩阵 间接读取

似乎不用特别的节省空间，尝试写一个浪费空间换取时间的函数

第三次尝试：取消rebuildcoord矩阵 直接读取



### 0731：

读取矩阵花费的时间比单独计算还要花时间 似乎

加速计算 还是 绕开计算

得到了一个可能在k增大的情况下，会有所加速的Sumdistance函数，同时代码变得简洁，适合做并行化。

多核异构并行计算 这本pdf写的很好 很入门

做完对sumdistance函数的算法优化之后，进行openmp并行优化：
对单个 chebyshevdistance  openmp  reduction max

### 0801:

尝试使用超算跑一次程序 开并行之后就始终段错误

初步定位是对于P2PDistance矩阵的malloc的问题,似乎这个东西不能开太大,虽然在并行线程里我只对这个矩阵进行了读操作

当n=500时,这个东西就已经很大了,还是说最开始的优化方向是错误的,还是分时间进行计算更加符合并行的优化

除此之外找了一下如何在超算上面开oneapi以及一众intel的编译器

source setvars.sh 脚本即可

开会总结:

- [x] 把icpc g++的问题报告
- [x] 修一修Sumdistance的bug

### 0804

check前几天不在的进度,使用gdb在服务器上跑了一遍,发现是因为最后的fclose函数的输入有问题,不知道为什么在本地跑程序的时候没有发现这个问题...但是最后还是找到了..

## 对SumDistance做并行优化, 直接优化:11300ms 删除rebuild版本:19000.ms

## 跑最大32线程 得到耗时:5000-6000ms左右

## 用ICPC编译 3700ms左右

### 0805

对Combination改循环

### 0806

有了一个初步的for循环遍历,由于数据依赖的关系,只能对大块的循环进行并行

小demo跑Case2的遍历要250000ms左右,先对这个进行大块的并行优化试一试

### 0807

改循环成功,还没有加并行,改成循环的耗时跟递归差不多

改并行的难点在于,怎么分配并行任务
$$
C^5_{100} = ∑^{99}_{i=4}C^4_{i}
$$

$$
C^n_{m} = ∑^{m-1}_{i=n}C^n_{i}
$$

现在的问题在于,对于每一个确定的集合点集,如何快速求出它对应的下标.

设[0,1,2,3,4]的下标为1 逐个递增,逢100进1

参赛的账号最多2节点128核

一节点 64核

单纯跑一遍遍历 有以下数据:
$$
Case1: C^2_{500} = 124,750|||0.654ms
$$

$$
Case2: C^5_{100} = 75,287,520|||273.539ms
$$

可见Case1 Case2 的计算倍数差不多400倍,于是我找了一个差不多也是400倍的数据:
$$
C^4_{1000} = 41,417,124,750|||120000ms
$$
收尾的一个小问题:

线程数量与计算次数之间的关系

合理分配尤为重要

### 0808

现在对整个题目进行总结:

得分点,或者说是难点:

1. 并行优化对确定支撑点集的切比雪夫距离的计算

   ...

1. 循环优化取得每一个确定的支撑点集的过程

   ....

2. 并行优化如何解决TopK问题

   ...

努力把第三点完成吧....

目前对第三点的思考,

堆排序+合并

时间复杂度应该算是	O(M+(n-M)logM)/core

最后的合并就可以

嵌套循环 32 2?

### 0809

开始写堆排

就三个函数,make_heap sort_heap   和 push_heap

写完了

case1快的不多,快了100ms左右 现在大概是3500-3600ms左右

ipcp 编译指令 -ipo过程优化 --fast-math浮点数计算优化  3200ms左右

绑了下核,3100ms左右

明天开始准备给Combination做并行优化

一直在思考的问题:并行netest之间的线程如何分配

### 0810

睡前跑32 64 case2

case2 32 2170842     700,000ms

case2 64 2170854	 900,000ms

case2 16 2170890	100,000ms

预测一下应该是64更加慢一点



的确是这样

对于n^2的循环得找到一个合适的线程数去匹配

少了一行数据....

case2 refer line 251

16 35 77 80 99 输出的result.txt没有这一行

### 0811

开始堆combination做并行

优先考虑的是大小嵌套循环怎么样才能得到一个比较好的效果

先写一个串行的test文件试一试吧

嵌套并行很慢很慢

所以放弃对sumdistance循环的并行,只对外部的combination做并行

```powershell
cu:

#查看CPU信息（型号）
cat /proc/cpuinfo | grep name | cut -f2 -d: | uniq -c
 32  Intel(R) Xeon(R) Gold 6326 CPU @ 2.90GHz
# 查看物理CPU个数
cat /proc/cpuinfo | grep "physical id" | sort | uniq | wc -l
2
# 查看每个物理CPU中core的个数(即核数)
cat /proc/cpuinfo | grep "cpu cores" | uniq
cpu cores: 16

# 查看逻辑CPU的个数
cat /proc/cpuinfo | grep "processor" | wc -l
32
```

写完了

开始测试

段错误 pivots

重新理一理思路 代码写的有点乱,有点烂

1.先从Combination开始,把flag和cnt优化一下

优化完成 ,用一个指针变量去做这个事

2.继续测试,意识到指针变量不能作为一个私有变量去使用,所以必须在并行块里面作malloc.但是,对于只读的P2PDist矩阵就可以直接作为一个share变量去使用.

```
icc pivot_0812.c -Ofast -ipo -qopenmp -fp-model fast=2 -g -o pivot_0812
```

3.继续测试,发现一开始有一个变量忘记改了,

pivot[i]=i

导致Startnum函数的输出矩阵有问题



vim翻半页

- `ctr-d`：向后翻半页
- `ctr-u`：向前翻半页

 vim整整页

- `ctr+f`：向后翻整页
- `ctr+b`：向前翻整页



有一个莫名其妙的bug

读不了共享变量StartPivots本地跑的可以读　cu上的就不行...

```shell
//Default version
gcc pivot.c -lm -o pivot
gcc pivot.c -lm -O3 -o pivot
//using intel icc compiler
icc pivot.c -Ofast -o pivot
//Add -fp-model to speedup floats caculation
//icc pivot.c -Ofast -fp-model fast=2 -o pivot
//Add -qopenmp to use openmp parallel
icc pivot.c -Ofast -fp-model fast=2 -qopenmp -o pivot
//Add ipo to reduce process costs
icc pivot.c -Ofast -fp-model fast=2 -qopenmp -ipo -o pivot
//bind proc to reduce cache conflict
export OMP_PROC_BIND=true
export OMP_PLACES=cores
```

继续昨天晚上的debug

使用gdb调试,发现值传递没有进行 很奇怪.

值传递的语句没有进行

```c
#define M 20000
#define N 5000

int main(int argc, char* argv[]){
    int i;
    int* B = (int*)malloc(sizeof(int)*M);
    for(i=0;i<M;i++) B[i]=i;
    #pragma omp parallel num_threads(4)
    {
        int tid = omp_get_thread_num();
        int* A = (int*)malloc(sizeof(int)*N);
        #pragma omp critical
        {
            for(i=0;i<N;i++){
                A[i] = tid*N+i;
            }
        }
        for(int i=0;i<N;i++){
			B[tid*N+i] = A[i];
        }

    }
    for(i=0;i<M;i++){
        printf("%d\n",B[i]);
    }

    free(B);
    return 0;
}
//输出正确 如果消去第一个的critical将会出错
```

### 0813

从头到尾捋了一遍代码

在Caculate_sort函数内出现了数据竞争,早上起来看看能不能解决

