### 出现的bug:

此pivot版本将先前的插入排序改为堆排序

跑case1正常通过,与result.txt文件比较后发现正常

跑case2的时候在maxdistancesum的环节出现了问题

问题点:

result.txt中缺少了refer-4dim-1h.txt中的line 251的一条数据

在排序前进行读取,发现计算部分没有问题,问题应该出现在之后的排序部分

复现问题的话,直接copy整个文件夹,在服务器上sbatch run.sh即可

