#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "util.h"

using namespace std;
using namespace Eigen;


class LPLC {
    public:
        LPLC(const vector<Vector3d>& vertices,const Matrix3Xi& faces, const vector<p_pair>& constraints);
        vector<Vector3d> deformed_; //变形后顶点
        void solveLinearSystem();
        void solveNonLinearSystem(int iterations);
        int iterations_ = 10; //非线性迭代次数
    private:
        vector<Vector3d> vertices_; //顶点
        vector<vector<int>> neighbors_; //邻居
        vector<p_pair> constraints_; //约束点对
        SparseMatrix<double> L_; //拉普拉斯矩阵
        vector<double> gamma_;  // |δ_i|
        vector<Vector3d> delta_hat_; // δ_i
        Eigen::Matrix3Xi g_faces;
        void buildLaplace();
        void buildNeighbors (const Matrix3Xi& faces);
        // non-liner
        void computeReferenceGamma();
        void computeDeltaHat();
        void solveWithRHS(const vector<Vector3d>& rhs);
};
