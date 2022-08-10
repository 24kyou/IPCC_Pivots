//
// Created by Admin on 2022/8/4.
//

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
float MAX(float a,float b){
    return a>b?a:b;
}
int main(){

//    part1 read file
using namespace  std;

    std::ifstream fin;
    int row=0,column=0,pivots=0;
     fin.open("./uniformvector-2dim-5h.txt");
     if(!fin.is_open()){
        fin.close();
        return 0;
     }
     fin>>column>>row>>pivots;
     Eigen::MatrixXd MatrixOrigin(row,column);//创建用于存储数据的矩阵
     for(int i =0;i<column;i++){
         for (int j=0;j<row;j++){
             fin>>MatrixOrigin(i,j);
         }
     }
    //计算每个点之间的距离 形成距离矩阵
    Eigen::MatrixXd MatrixDistance(row,row);//创建用于计算距离的矩阵
    Eigen::MatrixXd Maxtrix_new(1,column);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < row; ++j) {
//            MatrixDistance(i,j)=
            Maxtrix_new=MatrixOrigin.row(i)-MatrixOrigin.row(j);
            Maxtrix_new=Maxtrix_new.array().square();
            MatrixDistance(i,j)=sqrt(Maxtrix_new.sum());
        }
    }
    //释放不需要的矩阵
    MatrixOrigin.resize(0,0);
    Eigen::MatrixXd MatrixCheby(row,row);
    Eigen::MatrixXd MatrixCheby_new(row,row);
    for (int n = 0; n < pivots-1; ++n) {
//        打个标记在这里
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < row; ++j) {
                if(n==0){
                    Maxtrix_new=MatrixDistance.row(i)-MatrixDistance.row(j);
                    Maxtrix_new=Maxtrix_new.array().abs();

                }
            }
        }
    }
     fin.close();




}

