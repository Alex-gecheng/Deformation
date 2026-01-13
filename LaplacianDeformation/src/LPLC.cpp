#include <Eigen/Sparse>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include "LPLC.h"

using namespace std;
using namespace Eigen;

LPLC::LPLC(const vector<Vector3d>& vertices,  const Matrix3Xi& faces,const vector<p_pair>&  constraints)
    :vertices_(vertices),  constraints_(constraints) {
    deformed_ = vertices_; // 初始变形顶点为原始顶点
    g_faces = faces;
    buildNeighbors(faces);
    buildLaplace();
}

double cotan(const Vector3d &u, const Vector3d &v) {
    double denom = u.cross(v).norm();
    if (denom < 1e-12) return 0.0;
    return u.dot(v) / denom;
}
void LPLC::buildLaplace()
{
    typedef Triplet<double> T;
    int n = static_cast<int>(vertices_.size());
    vector<unordered_map<int,double>> weights(n);

    for (int f=0 ; f<g_faces.cols(); ++f) {
        int i = g_faces(0,f);
        int j = g_faces(1,f);
        int k = g_faces(2,f);

        const Vector3d &vi = vertices_[i];
        const Vector3d &vj = vertices_[j];
        const Vector3d &vk = vertices_[k];

        double cotA = cotan(vi-vj,vk-vi);
        double cotB = cotan(vj-vi,vk-vj);
        double cotC = cotan(vk-vj,vi-vk);

        weights[i][j] += 0.5 * cotC;
        weights[j][i] += 0.5 * cotC;
        weights[j][k] += 0.5 * cotA;
        weights[k][j] += 0.5 * cotA;
        weights[i][k] += 0.5 * cotB;
        weights[k][i] += 0.5 * cotB;
    }

    vector<T> triplets;
    for(int i = 0 ;i<n;i++){
        double diag =0.0;
        for(const auto &kv: weights[i]){
            int j = kv.first;
            double w = kv.second;
            triplets.emplace_back(i,j,-w);
            diag += w;
        }
        triplets.emplace_back(i,i,diag);
    }
    L_.resize(n, n);
    L_.setFromTriplets(triplets.begin(), triplets.end());
}

void LPLC::buildNeighbors (const Matrix3Xi& faces) {
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

void LPLC::solveLinearSystem() {
    int n = static_cast<int>(vertices_.size());
    if (n == 0) return;

    unordered_set<int> cset;
    unordered_map<int, Vector3d> cpos;
    for (const auto &c : constraints_) {
        cset.insert(c.first);
        cpos[c.first] = c.second;
    }

    if (cset.empty()) return; // nothing to anchor, system would be singular

    typedef Triplet<double> T;
    vector<T> triplets;
    triplets.reserve(L_.nonZeros());
    // 原始拉普拉斯坐标 delta = L * v_orig
    VectorXd ox(n), oy(n), oz(n);
    for (int i = 0; i < n; ++i) {
        ox[i] = vertices_[i].x();
        oy[i] = vertices_[i].y();
        oz[i] = vertices_[i].z();
    }

    VectorXd bx = L_ * ox;
    VectorXd by = L_ * oy;
    VectorXd bz = L_ * oz;

    // 约束点：覆盖对应行，改为 x_i = p_i
    for (int idx : cset) {
        triplets.emplace_back(idx, idx, 1.0);
        const Vector3d &p = cpos[idx];
        bx[idx] = p.x();
        by[idx] = p.y();
        bz[idx] = p.z();
    }

    // build system while moving constrained contributions to RHS
    for (int k = 0; k < L_.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(L_, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            double v = it.value();

            bool ci = cset.count(i);
            bool cj = cset.count(j);

            if (ci || cj) {
                if (!ci && cj) {
                    // 把 L_{ij} * p_j 移到右端项
                    const Vector3d &p = cpos[j];
                    bx[i] -= v * p.x();
                    by[i] -= v * p.y();
                    bz[i] -= v * p.z();
                }
                // skip adding matrix terms touching constraints (handled above)
                continue;
            }

            triplets.emplace_back(i, j, v);
        }
    }

    SparseMatrix<double> A(n, n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Success) {
        std::cerr << "Failed to factorize Laplacian system" << std::endl;
        return;
    }

    VectorXd x = solver.solve(bx);
    VectorXd y = solver.solve(by);
    VectorXd z = solver.solve(bz);
    if (solver.info() != Success) {
        std::cerr << "Failed to solve Laplacian system" << std::endl;
        return;
    }

    for (int i = 0; i < n; ++i) {
        deformed_[i] = Vector3d(x[i], y[i], z[i]);
    }
}



// non-liner
void LPLC::computeReferenceGamma()
{
    int n = vertices_.size();
    gamma_.resize(n);

    // 用面来累积每个顶点的一圈“扇形”向量和
    vector<Vector3d> accum(n, Vector3d::Zero());

    for (int f = 0; f < g_faces.cols(); ++f) {
        int i = g_faces(0, f);
        int j = g_faces(1, f);
        int k = g_faces(2, f);

        const Vector3d &vi = vertices_[i];
        const Vector3d &vj = vertices_[j];
        const Vector3d &vk = vertices_[k];

        // 以每个顶点为中心，累积相邻两条边的叉积
        accum[i] += (vj - vi).cross(vk - vi);
        accum[j] += (vk - vj).cross(vi - vj);
        accum[k] += (vi - vk).cross(vj - vk);
    }

    for (int i = 0; i < n; ++i) {
        gamma_[i] = accum[i].norm();
    }
}

void LPLC::computeDeltaHat()
{
    int n = vertices_.size();
    delta_hat_.assign(n, Vector3d::Zero());

    // 同样通过面来累积当前变形后的扇形向量
    vector<Vector3d> accum(n, Vector3d::Zero());

    for (int f = 0; f < g_faces.cols(); ++f) {
        int i = g_faces(0, f);
        int j = g_faces(1, f);
        int k = g_faces(2, f);

        const Vector3d &vi = deformed_[i];
        const Vector3d &vj = deformed_[j];
        const Vector3d &vk = deformed_[k];

        accum[i] += (vj - vi).cross(vk - vi);
        accum[j] += (vk - vj).cross(vi - vj);
        accum[k] += (vi - vk).cross(vj - vk);
    }

    for (int i = 0; i < n; ++i) {
        const Vector3d &d = accum[i];
        double norm_d = d.norm();
        if (norm_d > 1e-8 && gamma_[i] > 0.0) {
            delta_hat_[i] = gamma_[i] * d / norm_d;
        }
    }
}


void LPLC::solveWithRHS(const vector<Vector3d>& rhs){
    int n = static_cast<int>(vertices_.size());
    if (n == 0) return;

    unordered_set<int> cset;
    unordered_map<int, Vector3d> cpos;
    for (const auto &c : constraints_) {
        cset.insert(c.first);
        cpos[c.first] = c.second;
    }

    if (cset.empty()) return; // nothing to anchor, system would be singular

    typedef Triplet<double> T;
    vector<T> triplets;
    triplets.reserve(L_.nonZeros());
    // 原始拉普拉斯坐标 delta = L * v_orig
    VectorXd ox(n), oy(n), oz(n);
    VectorXd bx(n), by(n), bz(n);
    for (int i = 0; i < n; ++i) {
        bx[i] = rhs[i].x();
        by[i] = rhs[i].y();
        bz[i] = rhs[i].z();
    }

    // 约束点：覆盖对应行，改为 x_i = p_i
    for (int idx : cset) {
        triplets.emplace_back(idx, idx, 1.0);
        const Vector3d &p = cpos[idx];
        bx[idx] = p.x();
        by[idx] = p.y();
        bz[idx] = p.z();
    }

    // build system while moving constrained contributions to RHS
    for (int k = 0; k < L_.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(L_, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            double v = it.value();

            bool ci = cset.count(i);
            bool cj = cset.count(j);

            if (ci || cj) {
                if (!ci && cj) {
                    // 把 L_{ij} * p_j 移到右端项
                    const Vector3d &p = cpos[j];
                    bx[i] -= v * p.x();
                    by[i] -= v * p.y();
                    bz[i] -= v * p.z();
                }
                // skip adding matrix terms touching constraints (handled above)
                continue;
            }

            triplets.emplace_back(i, j, v);
        }
    }

    SparseMatrix<double> A(n, n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Success) {
        std::cerr << "Failed to factorize non-linear Laplacian system" << std::endl;
        return;
    }

    VectorXd x = solver.solve(bx);
    VectorXd y = solver.solve(by);
    VectorXd z = solver.solve(bz);
    if (solver.info() != Success) {
        std::cerr << "Failed to solve non-linear Laplacian system" << std::endl;
        return;
    }

    for (int i = 0; i < n; ++i) {
        deformed_[i] = Vector3d(x[i], y[i], z[i]);
    }
}

void LPLC::solveNonLinearSystem(int iterations)
{
    computeReferenceGamma();

    for (int it = 0; it < iterations; ++it) {
        computeDeltaHat();
        solveWithRHS(delta_hat_); 
    }
}
