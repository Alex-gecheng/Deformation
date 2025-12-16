#include "ARAP.h"
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <iostream>
#include <unordered_set>

using namespace std;
using namespace Eigen;

ARAP::ARAP(const vector<Vector3d>& vertices,  const Matrix3Xi& faces,const vector<p_pair>&  constraints)
    :vertices_(vertices),  constraints_(constraints) {
    deformed_ = vertices_; // 初始变形顶点为原始顶点
    buildNeighbors(faces);
    rotations_.resize(vertices.size(), Matrix3d::Identity()); // 初始化旋转矩阵为单位矩阵
    buildLaplace();
}


void  ARAP::run(int iterations) {
    for(int i=0; i<iterations; i++){
        LocalStep();
        GlobalStep();
    }
}

void ARAP::buildLaplace()
{
    typedef Triplet<double> T;
    vector<T> triplets;

    int n = vertices_.size();
    for (int i = 0; i < n; ++i) {
        double wsum = 0.0;
        for (int j : neighbors_[i]) {
            triplets.push_back(T(i, j, -1.0)); // 均一权重
            wsum += 1.0;
        }
        triplets.push_back(T(i, i, wsum)); // 对角元素
    }

    L_.resize(n, n);
    L_.setFromTriplets(triplets.begin(), triplets.end());
}

void ARAP::LocalStep()
{
    for (size_t i = 0; i < vertices_.size(); ++i) {
        Matrix3d S = Matrix3d::Zero();

        for (int j : neighbors_[i]) {
            Vector3d pi = vertices_[i];
            Vector3d pj = vertices_[j];
            Vector3d qi = deformed_[i];
            Vector3d qj = deformed_[j];

            // 计算局部协方差矩阵
            S += (pj - pi) * (qj - qi).transpose();
        }

        // SVD 分解
        JacobiSVD<Matrix3d> svd(S, ComputeFullU | ComputeFullV);
        Matrix3d R = svd.matrixV() * svd.matrixU().transpose();

        // 确保旋转矩阵是正交的
        if (R.determinant() < 0) {
            Matrix3d V = svd.matrixV();
            V.col(2) *= -1;
            R = V * svd.matrixU().transpose();
        }

        rotations_[i] = R;
    }
}


void ARAP::GlobalStep()
{
    const int n = vertices_.size();

    VectorXd bx = VectorXd::Zero(n);
    VectorXd by = VectorXd::Zero(n);
    VectorXd bz = VectorXd::Zero(n);

    /* ===============================
       1. 构建 ARAP 右端项 b
       b_i = sum_j 0.5 * (R_i + R_j) * (p_i - p_j)
       =============================== */
    for (int i = 0; i < n; ++i) {
        for (int j : neighbors_[i]) {
            Vector3d pij = vertices_[i] - vertices_[j];
            Vector3d term = 0.5 * (rotations_[i] + rotations_[j]) * pij;

            bx(i) += term.x();
            by(i) += term.y();
            bz(i) += term.z();
        }
    }

    //2. 复制 Laplacian，准备加硬约束
    SparseMatrix<double> A = L_;

    /* ===============================
       3. Dirichlet 约束（硬约束）
       对每个约束点 i:
       - 清空 A 的第 i 行
       - A(i,i) = 1
       - b(i) = target
       =============================== */
    for (const auto& c : constraints_) {
    int i = c.first;
    const Vector3d& target = c.second;

    assert(i >= 0 && i < n);

    // 清空第 i 行（注意：SparseMatrix 默认是列主序）
    for (int col = 0; col < A.outerSize(); ++col) {
        for (SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            if (it.row() == i) {
                it.valueRef() = 0.0;
            }
        }
    }

    A.coeffRef(i, i) = 1.0;

    bx(i) = target.x();
    by(i) = target.y();
    bz(i) = target.z();
}

    //4. 求解三个线性系统
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(A);

    if (solver.info() != Success) {
        std::cerr << "ARAP GlobalStep: decomposition failed!" << std::endl;
        return;
    }

    VectorXd new_x = solver.solve(bx);
    VectorXd new_y = solver.solve(by);
    VectorXd new_z = solver.solve(bz);

    if (solver.info() != Success) {
        std::cerr << "ARAP GlobalStep: solve failed!" << std::endl;
        return;
    }

    //5. 更新变形顶点
    for (int i = 0; i < n; ++i) {
        deformed_[i].x() = new_x(i);
        deformed_[i].y() = new_y(i);
        deformed_[i].z() = new_z(i);
    }
}

void ARAP::buildNeighbors (const Matrix3Xi& faces) {
    int n = vertices_.size();
    vector<unordered_set<int>> temp_neighbors(n);

    for (int f = 0; f < faces.cols(); ++f) {
        int i = faces(0,f);
        int j = faces(1,f);
        int k = faces(2,f);

        temp_neighbors[i].insert(j); temp_neighbors[i].insert(k);
        temp_neighbors[j].insert(i); temp_neighbors[j].insert(k);
        temp_neighbors[k].insert(i); temp_neighbors[k].insert(j);
    }

    neighbors_.resize(n);
    for (int v = 0; v < n; ++v) {
        neighbors_[v] = vector<int>(temp_neighbors[v].begin(), temp_neighbors[v].end());
    }
}
