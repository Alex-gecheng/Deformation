#include <Eigen/Sparse>
#include <iostream>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include "LPLC.h"

using namespace std;
using namespace Eigen;

LPLC::LPLC(const vector<Vector3d>& vertices,  const Matrix3Xi& faces,const vector<p_pair>&  constraints, SparseMatrix<double> Gamma)
    :vertices_(vertices),  constraints_(constraints), Gamma_(Gamma) {
    deformed_ = vertices_; // 初始变形顶点为原始顶点
    faces_ = faces;
    buildNeighbors();
    buildLaplace();
    computeAi();
    computemu();
}

boneSegment::boneSegment(const vector<int>& controller_indices, const vector<Vector3d>& ab_points, const vector<Vector3d>& vertices,int r){
    user_controller = controller_indices;
    user_ab = ab_points;
    r_ = r;
    vertices_ = vertices;
    buildBox();
    computeKandTheta();
    buildGamma();
}

void LPLC::buildLaplace(){
    // 找对角点
    map <pair<int,int>, vector<int>> edge_faces;
    for (int i = 0 ; i < faces_.cols(); ++i){
        int v0 = faces_(0,i);
        int v1 = faces_(1,i);
        int v2 = faces_(2,i);

        edge_faces[{min(v0,v1),max(v0,v1)}].push_back(v2);
        edge_faces[{min(v1,v2),max(v1,v2)}].push_back(v0);
        edge_faces[{min(v2,v0),max(v2,v0)}].push_back(v1);
    }

    // 遍历邻居，计算w    i , <j,w>
    weights_.resize(vertices_.size());

    for (int j = 0 ; j < vertices_.size(); ++j){
        for (int k = 0 ; k < neighbors_[j].size(); ++k){
            int neighbor = neighbors_[j][k];
            for (auto &opposite : edge_faces[{min(j,neighbor),max(j,neighbor)}]){
                // 两条邻边
                Vector3d v1 = vertices_[opposite] - vertices_[neighbor];
                Vector3d v2 = vertices_[opposite] - vertices_[j];
                double cot_angle = v1.dot(v2) / (v1.cross(v2)).norm();
                // 存储权重
                weights_[j][neighbor] += cot_angle;
            }
            weights_[j][neighbor] *= 0.5;
            weights_[neighbor][j] = weights_[j][neighbor]; // 对称
        }
    }

    // 构建L_
    L_.resize(vertices_.size(), vertices_.size());
    for (int i = 0; i < vertices_.size(); ++i){
        for (auto &kv : weights_[i]){
            int j = kv.first;
            double w = kv.second;
            L_.insert(i,j) = -w;
        }
        double diag = 0.0;
        for (auto &kv : weights_[i]){
            diag += kv.second;
        }
        L_.insert(i,i) = diag;
    }
    L_.makeCompressed();

}

void LPLC::buildNeighbors(){
    neighbors_.resize(vertices_.size());
    for (int i = 0 ; i<faces_.cols();++ i){
        int v0 = faces_(0,i);
        int v1 = faces_(1,i);
        int v2 = faces_(2,i);

        neighbors_[v0].push_back(v1);
        neighbors_[v0].push_back(v2);
        neighbors_[v1].push_back(v0);
        neighbors_[v1].push_back(v2);
        neighbors_[v2].push_back(v0);
        neighbors_[v2].push_back(v1);
    }

    for(auto &one : neighbors_){
        sort(one.begin(),one.end());
        one.erase(unique(one.begin(),one.end()),one.end());
    }
    
}

void LPLC::computeAi() {
    int n = vertices_.size();
    Ai_.resize(n); // 每个顶点存 Ai 矩阵，列数 = 邻居数

    for (int i = 0; i < n; ++i) {
        int ni = neighbors_[i].size(); 
        if (ni < 3) continue; 
        Matrix3Xd Ai(3, ni);

        for (int j = 0; j < ni; ++j) {
            int prev = neighbors_[i][(j + ni - 1) % ni];
            int cur  = neighbors_[i][j];
            Vector3d normal = (vertices_[prev] - vertices_[i]).cross(vertices_[cur] - vertices_[i]);
            Ai.col(j) = normal;
        }
        Ai_[i] = Ai;
    }
}

void LPLC::computemu() {
    int n = vertices_.size();
    mu_.resize(n); 

    vector<Vector3d> delta(n, Vector3d::Zero());
    for (int i = 0; i < n; ++i) {
        for (int j : neighbors_[i]) {
            double w = weights_[i][j];
            delta[i] += w * (vertices_[i] - vertices_[j]);
        }
    }
    delta_ = delta; 
    for (int i = 0; i < n; ++i) {
        // 跳过边界点或 Ai 未构建的点
        if (Ai_[i].cols() < 3) {
            mu_[i] = VectorXd::Zero(Ai_[i].cols());
            continue;
        }

        Matrix3Xd& Ai = Ai_[i]; 
        Vector3d delta_i = delta[i];

        // 计算 μ_i * Ai= δ_i
        JacobiSVD<Matrix3Xd> svd(Ai, ComputeThinU | ComputeThinV);
        VectorXd mu_i = svd.solve(delta_i);
        mu_[i] = mu_i;
    }
}


void LPLC::computeDix(){
    int n = vertices_.size();
    Dix_.resize(n);
    for (int i = 0; i < n; ++i) {
        int ni = neighbors_[i].size();
        if (ni < 3) {
            continue;
        }
        Vector3d Di;
        Di = Vector3d::Zero();
        for (int j = 0; j < ni; ++j) {
            int prev = neighbors_[i][(j + ni - 1) % ni];
            int cur  = neighbors_[i][j];

            Vector3d normal = (deformed_[prev] - deformed_[i]).cross(deformed_[cur] - deformed_[i]);

            Di += mu_[i][j] * normal;
        }
        const double eps = 1e-8;
        if (Di.norm() < eps)
            Dix_[i] = Vector3d::Zero();
        else
            Dix_[i] = Di * delta_[i].norm() / Di.norm();
    }
}


// 约束盒子    用户给四个点，获取范围内所有点
void boneSegment::buildBox(){
    vector<Vector3d> user_controller_vertices;
    for (int i = 0; i<user_controller.size(); ++i){
        user_controller_vertices.push_back( vertices_[ user_controller[i] ] );
    }
    double x_min = user_controller_vertices[0].x();
    double x_max = user_controller_vertices[0].x();
    double y_min = user_controller_vertices[0].y();
    double y_max = user_controller_vertices[0].y();
    double z_min = user_controller_vertices[0].z();
    double z_max = user_controller_vertices[0].z();
    for(auto &v : user_controller_vertices){
        if (v.x() < x_min) x_min = v.x();
        if (v.x() > x_max) x_max = v.x();
        if (v.y() < y_min) y_min = v.y();
        if (v.y() > y_max) y_max = v.y();
        if (v.z() < z_min) z_min = v.z();
        if (v.z() > z_max) z_max = v.z();
    }
    for (int j = 0; j < vertices_.size(); ++j){
        const Vector3d& v = vertices_[j];
        if (v.x() >= x_min - r_ && v.x() <= x_max + r_ &&
            v.y() >= y_min - r_ && v.y() <= y_max + r_ &&
            v.z() >= z_min - r_ && v.z() <= z_max + r_ ){
            box_.push_back(j);  
        }
    }

    // 构建采样点骨架
    Vector3d ab = user_ab[1] - user_ab[0];
    limit_vertices_.reserve(r_);
    for (int i = 0; i < r_; ++i){
        double t = double(i) / double(r_ - 1);
        Vector3d sample_point = user_ab[0] + t * ab;
        limit_vertices_.push_back(sample_point);
    }
}

void boneSegment::computeKandTheta(){
    K_.resize(limit_vertices_.size());
    for (int i = 0; i < limit_vertices_.size(); ++i){
        Vector3d Si = limit_vertices_[i];
        double sum = 0.0;
        for (int j : box_){
            Vector3d diff = Si - vertices_[j];
            double dist2 = diff.squaredNorm() + 1e-8; 
            double weight = 1.0 / dist2;
            K_[i][j]= weight;
            sum += weight;
        }
        for (auto& kv : K_[i]){
            kv.second /= sum;
        }
    }
}

void boneSegment::buildGamma(){
    Gamma_.resize(r_-1, vertices_.size());
    for (int i = 1; i < r_; ++i){
        for (int j : box_){
            double k_ij = K_[i][j];
            double k_ij1 = K_[i-1][j];
            double k_r = K_[r_-1][j];
            double k_0 = K_[0][j];
            double val = (k_ij - k_ij1) - 1.0/(r_-1)*(k_r - k_0);
            Gamma_.coeffRef(i-1,j) = val;
        }
    }
     Gamma_.makeCompressed(); 
}


// void LPLC::solve(int iterations)
// {
//     int n = vertices_.size();
//     deformed_ = vertices_; // 初始化
//     int m = Gamma_.rows();
//     for (int i = 0; i < iterations; ++i)
//     {
//         computeDix();
//         VectorXd bx(n), by(n), bz(n);
//         for (int i = 0; i < n; ++i) {
//             bx[i] = Dix_[i].x();
//             by[i] = Dix_[i].y();
//             bz[i] = Dix_[i].z();
//         }

//         //constraints
//         SparseMatrix<double> A = L_;
//         VectorXd bx_c = bx, by_c = by, bz_c = bz;

//         // 约束  LV = delta   修改L
//         // for (const auto &p : constraints_){
//         //     int i = p.first;

//         //     for (int j = 0; j < vertices_.size(); ++j){
//         //         if (j != i) L_.coeffRef(i,j) = 0.0;
//         //     }
//         //     L_.coeffRef(i,i) = 1.0;

//         //     bx_c[i] = p.second.x();
//         //     by_c[i] = p.second.y();
//         //     bz_c[i] = p.second.z();
//         // }

//         // SparseLU<SparseMatrix<double>> solver;
//         // solver.compute(A);

//         // VectorXd Xx = solver.solve(bx_c);
//         // VectorXd Xy = solver.solve(by_c);
//         // VectorXd Xz = solver.solve(bz_c);


//         int rows_aug = A.rows() + m;
//         SparseMatrix<double> A_aug(rows_aug, n);
//         vector<Triplet<double>> trips;

//         // L 部分
//         for (int k = 0; k < A.outerSize(); ++k)
//             for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
//                 trips.emplace_back(it.row(), it.col(), it.value());

//         // Γ 部分
//         double s = std::sqrt(lambda_skel_);
//         for (int k = 0; k < Gamma_.outerSize(); ++k)
//             for (SparseMatrix<double>::InnerIterator it(Gamma_, k); it; ++it)
//                 trips.emplace_back(A.rows() + it.row(), it.col(), s * it.value());

//         A_aug.setFromTriplets(trips.begin(), trips.end());

//         VectorXd bx_aug(rows_aug), by_aug(rows_aug), bz_aug(rows_aug);
//         bx_aug.setZero();
//         by_aug.setZero();
//         bz_aug.setZero();

//         // 原 Laplacian 右端
//         bx_aug.head(n) = bx_c;
//         by_aug.head(n) = by_c;
//         bz_aug.head(n) = bz_c;

//         for (const auto &p : constraints_){
//             int i = p.first;

//             for (int j = 0; j < vertices_.size(); ++j){
//                 if (j != i) L_.coeffRef(i,j) = 0.0;
//             }
//             L_.coeffRef(i,i) = 1.0;

//             bx_aug[i] = p.second.x();
//             by_aug[i] = p.second.y();
//             bz_aug[i] = p.second.z();
//         }

//         // 1. 构造正规方程
//         SparseMatrix<double> AtA = A_aug.transpose() * A_aug;

//         VectorXd Atbx = A_aug.transpose() * bx_aug;
//         VectorXd Atby = A_aug.transpose() * by_aug;
//         VectorXd Atbz = A_aug.transpose() * bz_aug;

//         // 2. 用 SparseLU 解方阵
//         SparseLU<SparseMatrix<double>> solver;
//         solver.compute(AtA);

//         if (solver.info() != Success) {
//             std::cerr << "Decomposition failed!" << std::endl;
//             return;
//         }

//         VectorXd Xx = solver.solve(Atbx);
//         VectorXd Xy = solver.solve(Atby);
//         VectorXd Xz = solver.solve(Atbz);
//         for (int j = 0; j < n; ++j){
//             deformed_[j] = Vector3d(Xx[j], Xy[j], Xz[j]);
//         }
            
//     }
// }

void LPLC::solve(int iterations)
{
    int n = vertices_.size();          // 顶点数
    int m = Gamma_.rows();             // 骨骼约束行数
    int c = constraints_.size();       // 固定点数量

    deformed_ = vertices_;

    for (int it = 0; it < iterations; ++it)
    {
        // 1. 计算 Laplacian 右端
        computeDix();

        VectorXd bx(n), by(n), bz(n);
        for (int i = 0; i < n; ++i) {
            bx[i] = Dix_[i].x();
            by[i] = Dix_[i].y();
            bz[i] = Dix_[i].z();
        }
        // 2. 构建固定点矩阵 C
        SparseMatrix<double> C(c, n);
        VectorXd cx(c), cy(c), cz(c);

        int row = 0;
        for (auto &p : constraints_) {
            int vi = p.first;
            C.insert(row, vi) = 1.0;
            cx[row] = p.second.x();
            cy[row] = p.second.y();
            cz[row] = p.second.z();
            row++;
        }

        // 3. 拼接 A_aug = [ L ; sqrt(lambda)*Γ ; C ]
        int rows_aug = n + m + c;
        SparseMatrix<double> A_aug(rows_aug, n);
        vector<Triplet<double>> trips;

        // --- L ---
        for (int k = 0; k < L_.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(L_, k); it; ++it)
                trips.emplace_back(it.row(), it.col(), it.value());

        // --- Γ ---
        double s = std::sqrt(lambda_skel_);
        int off1 = n;
        for (int k = 0; k < Gamma_.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(Gamma_, k); it; ++it)
                trips.emplace_back(off1 + it.row(), it.col(), s * it.value());

        // --- C ---
        int off2 = off1 + m;
        for (int k = 0; k < C.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(C, k); it; ++it)
                trips.emplace_back(off2 + it.row(), it.col(), it.value());

        A_aug.setFromTriplets(trips.begin(), trips.end());
        A_aug.makeCompressed();

        // 4. 拼右端
        VectorXd bx_aug(rows_aug), by_aug(rows_aug), bz_aug(rows_aug);
        bx_aug.setZero();
        by_aug.setZero();
        bz_aug.setZero();

        // Laplacian
        bx_aug.head(n) = bx;
        by_aug.head(n) = by;
        bz_aug.head(n) = bz;

        // Gamma 是 0
        // 固定点
        bx_aug.segment(off2, c) = cx;
        by_aug.segment(off2, c) = cy;
        bz_aug.segment(off2, c) = cz;

        // ========================================
        // 5. 正规方程
        // ========================================
        SparseMatrix<double> AtA = A_aug.transpose() * A_aug;

        VectorXd Atbx = A_aug.transpose() * bx_aug;
        VectorXd Atby = A_aug.transpose() * by_aug;
        VectorXd Atbz = A_aug.transpose() * bz_aug;

        // ========================================
        // 6. 解线性系统
        // ========================================
        SparseLU<SparseMatrix<double>> solver;
        solver.compute(AtA);

        if (solver.info() != Success) {
            std::cerr << "Decomposition failed!" << std::endl;
            return;
        }

        VectorXd Xx = solver.solve(Atbx);
        VectorXd Xy = solver.solve(Atby);
        VectorXd Xz = solver.solve(Atbz);

        // ========================================
        // 7. 更新顶点
        // ========================================
        for (int i = 0; i < n; ++i)
            deformed_[i] = Vector3d(Xx[i], Xy[i], Xz[i]);
    }
}