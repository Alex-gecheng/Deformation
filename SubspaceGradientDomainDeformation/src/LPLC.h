#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "util.h"

using namespace std;
using namespace Eigen;


class LPLC {
    public:
        LPLC(const vector<Vector3d>& vertices,const Matrix3Xi& faces, const vector<p_pair>& constraints ,SparseMatrix<double> Gamma);
        vector<Vector3d> deformed_; //变形后顶点
        void solve(int iterations);
        int iterations_ = 10; //非线性迭代次数
    private:
        vector<Vector3d> vertices_;
        SparseMatrix<double> Gamma_;
        Matrix3Xi faces_;
        vector<vector<int>> neighbors_; 
        vector<p_pair> constraints_; 
        SparseMatrix<double> L_; //拉普拉斯矩阵
        vector<map<int,double>> weights_; 
        vector<Vector3d> delta_;
        vector<Matrix3Xd> Ai_; 
        vector<VectorXd> mu_;
        vector<Vector3d> Dix_;
        double lambda_skel_ = 1000;
        void buildLaplace();
        void buildNeighbors ();
        void computeAi();
        void computemu();
        void computeDix();
};


class boneSegment{
            public:
                boneSegment(const vector<int>& controller_indices, const vector<Vector3d>& ab_points,const vector<Vector3d>& vertices,int r=9);
                vector<int> user_controller;   //36,41,191,186
                vector<Vector3d> user_ab; //用户直接给出两个端点  (-10,0,0) (-2,0,0)
                vector<int> box_;
                int r_ = 9;    //采样点数量
                vector<Vector3d> limit_vertices_; //采样点Si
                vector<map<int,double>> K_;   // K_[i][j] = k_ij  Si= KX K=omiga/sum
                vector<double> Theta_;
                SparseMatrix<double> Gamma_; 
            private:
                vector<Vector3d> vertices_;
                void buildBox();
                void computeKandTheta();
                void buildGamma();
};