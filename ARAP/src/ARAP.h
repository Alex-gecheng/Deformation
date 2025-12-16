#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "util.h"

using namespace std;
using namespace Eigen;


class ARAP {
    public:
        ARAP(const vector<Vector3d>& vertices,const Matrix3Xi& faces, const vector<p_pair>& constraints);
        void run(int iterations = 10);
        vector<Vector3d> deformed_; //变形后顶点
    private:
        vector<Vector3d> vertices_; //顶点
        vector<vector<int>> neighbors_; //邻居
        vector<p_pair> constraints_; //约束点对
        vector<Matrix3d> rotations_; //局部旋转矩阵
        SparseMatrix<double> L_; 
        void LocalStep();
        void GlobalStep();
        void buildLaplace();
        void buildNeighbors (const Matrix3Xi& faces);

};
