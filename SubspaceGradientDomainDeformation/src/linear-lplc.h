#pragma once 
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "util.h"

using namespace std;
using namespace Eigen;

class linearLPLC {
    public:
        linearLPLC(const vector<Vector3d>& vertices, const Matrix3Xi& faces, const vector<p_pair>& constraints);
        vector<Vector3d> deformed_; //变形后顶点
        void solveLinearSystem();
    private:
        vector<Vector3d> vertices_; //顶点
        Matrix3Xi faces_;
        vector<vector<int>> neighbors_; //邻居
        vector<p_pair> constraints_; //约束点对
        SparseMatrix<double> L_; //拉普拉斯矩阵
        vector<map<int,double>> weights_; //权重
        void buildNeighbors();
        void setConstraints();
        void buildLaplace();
};