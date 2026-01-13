#include "linear-lplc.h"

linearLPLC::linearLPLC(const vector<Vector3d>& vertices, const Matrix3Xi& faces, const vector<p_pair>& constraints) {
    vertices_ = vertices;
    constraints_ = constraints;
    deformed_ = vertices_; // 初始变形顶点为原始顶点
    faces_= faces;
    buildNeighbors();
    buildLaplace();
}


void linearLPLC::buildNeighbors(){
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

void linearLPLC::buildLaplace(){

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
                // 计算cot值
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
void linearLPLC::solveLinearSystem(){
    // 计算原始 Laplacian 坐标
    VectorXd delta_x(vertices_.size());
    VectorXd delta_y(vertices_.size());
    VectorXd delta_z(vertices_.size());
    delta_x.setZero();
    delta_y.setZero();
    delta_z.setZero();
    for (int i = 0 ; i < vertices_.size(); ++i){
        for (int j : neighbors_[i]){
            double w = weights_[i][j];
            delta_x[i] += w * (vertices_[i].x() - vertices_[j].x());
            delta_y[i] += w * (vertices_[i].y() - vertices_[j].y());
            delta_z[i] += w * (vertices_[i].z() - vertices_[j].z());
        }
    }
    // 约束  LV = delta   修改L
    for (const auto &p : constraints_){
        int i = p.first;

        for (int j = 0; j < vertices_.size(); ++j){
            if (j != i) L_.coeffRef(i,j) = 0.0;
        }
        L_.coeffRef(i,i) = 1.0;

        delta_x[i] = p.second.x();
        delta_y[i] = p.second.y();
        delta_z[i] = p.second.z();
    }

    SparseLU<SparseMatrix<double>> solver;
    solver.compute(L_);

    VectorXd Vx = solver.solve(delta_x);
    VectorXd Vy = solver.solve(delta_y);
    VectorXd Vz = solver.solve(delta_z);

    // 填充 deformed_ 顶点
    deformed_.resize(vertices_.size());
    for (int i = 0; i < vertices_.size(); ++i) {
        deformed_[i] = Vector3d(Vx[i], Vy[i], Vz[i]);
    }

}

