#include "mips.h"
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <iostream>

MIPS::MIPS(const Eigen::Matrix3Xd& Vin,
           const Eigen::Matrix3Xi& Fin,
           const std::vector<Constraint>& constraints_in)
    : V(Vin),
      F(Fin),
      nV(static_cast<int>(Vin.cols())),
      nF(static_cast<int>(Fin.cols())),
      constraints(constraints_in) {
    is_fixed.assign(nV, false);
    constraint_initial_uv.assign(nV, Eigen::Vector2d::Zero());
    constraint_active_uv.assign(nV, Eigen::Vector2d::Zero());
    constraint_target_uv.assign(nV, Eigen::Vector2d::Zero());

    for (int i = 0; i < nV; ++i) {
        Eigen::Vector2d p(V(0, i), V(1, i));
        constraint_initial_uv[i] = p;
        constraint_active_uv[i] = p;
        constraint_target_uv[i] = p;
    }

    for (const auto& c : constraints) {
        int vid = c.first;
        if (vid >= 0 && vid < nV) {
            is_fixed[vid] = true;
            constraint_target_uv[vid] = c.second.head<2>();
        }
    }

    buildLocal2DReference(2);
    initializeUVWithBoundary();
    buildVertexAdjacency();
    buildVertexFaceAdjacency();
    buildGraphColoring();
    total_energy = computeTotalEnergy();
}


bool MIPS::buildLocal2DReference(int dim) {
    face_data.resize(nF);
    if(dim != 2 && dim != 3) {
        std::cerr << "Unsupported dimension for local reference: " << dim << std::endl;
        return false;
    }
    if (dim == 2) {
        // 2D输入，直接使用顶点坐标作为局部参考
        for (int f = 0; f < nF; ++f) {
            int i = F(0, f);
            int j = F(1, f);
            int k = F(2, f);

            Eigen::Vector2d p0, p1, p2;    
            p0 = V.col(i).head<2>();     //每个顶点的x,y
            p1 = V.col(j).head<2>();
            p2 = V.col(k).head<2>();    

            Eigen::Matrix2d Dm;
            Dm.col(0) = p1 - p0;
            Dm.col(1) = p2 - p0;

            face_data[f].vid = Eigen::Vector3i(i, j, k);
            face_data[f].Dm = Dm;
            face_data[f].Dm_inv = Dm.inverse();
            face_data[f].area = 0.5 * std::abs(Dm.determinant());
            face_data[f].valid = true;
        }
        return true;
    }
    else{
        for (int f = 0; f < nF; ++f) {
            int i = F(0, f);
            int j = F(1, f);
            int k = F(2, f);

            Eigen::Vector3d p0, p1, p2;
            p0 = V.col(i);
            p1 = V.col(j);  
            p2 = V.col(k);

            Eigen::Vector3d e1 = p1 - p0;
            double len = e1.norm();
            if (len < 1e-12) {
                face_data[f].valid = false;
                continue;
            }
            e1 /= len;

            Eigen::Vector3d v2 = p2 - p0;
            double x2 = v2.dot(e1);
            Eigen::Vector3d ortho = v2 - x2 * e1;
            double y2 = ortho.norm();

            if (y2 < 1e-12) {
                face_data[f].valid = false;
                continue;
            }

            Eigen::Matrix2d Dm;
            Dm.col(0) = Eigen::Vector2d(len, 0.0);
            Dm.col(1) = Eigen::Vector2d(x2, y2);

            face_data[f].vid = Eigen::Vector3i(i, j, k);
            face_data[f].Dm = Dm;
            face_data[f].Dm_inv = Dm.inverse();
            face_data[f].area = 0.5 * std::abs(Dm.determinant());
            face_data[f].valid = true;
        }

        return true;
    }
    
}


bool MIPS::initializeUVWithBoundary() {
    UV.resize(2, nV);
    for (int i = 0; i < nV; ++i) {
        UV(0, i) = V(0, i);
        UV(1, i) = V(1, i);
    }
    applyConstraints();
    return true;
}

// 构建顶点邻接表
void MIPS::buildVertexAdjacency(){
    vertex_adj.resize(nV);
    for (int f = 0; f < nF; ++f) {
        int i = F(0,f);
        int j = F(1,f);
        int k = F(2,f);

        vertex_adj[i].push_back(j);
        vertex_adj[i].push_back(k);

        vertex_adj[j].push_back(i);
        vertex_adj[j].push_back(k);

        vertex_adj[k].push_back(i);
        vertex_adj[k].push_back(j);
    }
    for (int v = 0; v < nV; ++v) {
        auto& vec = vertex_adj[v];
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    }
}

// 构建面邻接表
void MIPS::buildVertexFaceAdjacency(){
    vertex_faces.resize(nV);
    for (int f = 0; f < nF; ++f) {
        int i = F(0,f);
        int j = F(1,f);
        int k = F(2,f);

        vertex_faces[i].push_back(f);
        vertex_faces[j].push_back(f);
        vertex_faces[k].push_back(f);
    }
}

// 图着色
void MIPS::buildGraphColoring() {
    std::vector<int> vertex_color(nV, -1);
    for (int v= 0 ; v < nV ;++ v){
        std::unordered_set<int> used_colors;
        for (auto adj : vertex_adj[v]) {
            if (vertex_color[adj] != -1) {
                used_colors.insert(vertex_color[adj]);
            }
            
        }
        int c = 0;
        while (used_colors.count(c)!=0) {
            c++;
        }
        vertex_color[v] = c;
        num_colors = std::max(num_colors, c+1);
    }

    // 分组
    int max_color = *std::max_element(vertex_color.begin(), vertex_color.end());

    color_groups.clear();
    color_groups.resize(num_colors);
    std::cout<<"color group num: "<<max_color+1<< std::endl;
    for (int v = 0; v < nV; ++v) {
        color_groups[vertex_color[v]].push_back(v);
    }
}

// 应用当前约束位置到 UV 中（每次更新 UV 后都要调用）
void MIPS::applyConstraints() {
    for (const auto& c : constraints) {
        int vid = c.first;
        if (vid >= 0 && vid < nV) {
            UV.col(vid) = constraint_active_uv[vid];
        }
    }
}

// 检查是否有翻转面存在
bool MIPS::hasAnyFlippedFace() const {
    for (int f = 0; f < nF; ++f) {
        if (isFaceFlipped(f)) {
            return true;
        }
    }
    return false;
}

// 约束点在初始值和目标值之间插值，分步推进约束生效
void MIPS::setConstraintProgress(double t) {
    t = std::max(0.0, std::min(1.0, t));
    for (const auto& c : constraints) {
        int vid = c.first;
        if (vid >= 0 && vid < nV) {
            constraint_active_uv[vid] =
                (1.0 - t) * constraint_initial_uv[vid] + t * constraint_target_uv[vid];
        }
    }
}

// 尝试推进约束到目标位置的某个阶段 t，若过程中出现翻转则自动回退步长，保证推进安全
bool MIPS::advanceConstraintProgressSafely(double target_t,
                                           int max_backtracks,
                                           double shrink) {
    target_t = std::max(constraint_progress, std::min(1.0, target_t));
    double cand_t = target_t;

    for (int k = 0; k < max_backtracks; ++k) {
        setConstraintProgress(cand_t);
        applyConstraints();
        if (!hasAnyFlippedFace()) {
            constraint_progress = cand_t;
            return true;
        }
        cand_t = constraint_progress + (cand_t - constraint_progress) * shrink;
    }

    setConstraintProgress(constraint_progress);
    applyConstraints();
    return false;
}

// 计算单个面片的雅可比矩阵 J = Ds * Dm^{-1}
Eigen::Matrix2d MIPS::computeJacobian(int fid) const {
    const auto& fd = face_data[fid];
    const auto& vid = fd.vid;

    Eigen::Vector2d u0 = UV.col(vid[0]);
    Eigen::Vector2d u1 = UV.col(vid[1]);
    Eigen::Vector2d u2 = UV.col(vid[2]);

    Eigen::Matrix2d Ds;
    Ds.col(0) = u1 - u0;
    Ds.col(1) = u2 - u0;

    return Ds * fd.Dm_inv;
}

// 计算单个面片的 MIPS 能量，E = ||J||^2 / det(J)，其中 J = Ds * Dm^{-1}
double MIPS::computeFaceEnergy(int fid) const {
    if (!face_data[fid].valid) {
        return 0.0;
    }

    Eigen::Matrix2d J = computeJacobian(fid);
    double detJ = J.determinant();
    if (detJ <= 1e-12) {
        return 1e20;
    }

    double frob2 = J.squaredNorm();
    return frob2 / detJ;
}

// 总能量
double MIPS::computeTotalEnergy() const {
    double total = 0.0;
    for (int f = 0; f < nF; ++f) {
        total += computeFaceEnergy(f)*face_data[f].area;    //乘以面积权重
    }
    return total;
}

// 相邻面能量
double MIPS::computeLocalEnergy(int v) const {
    double E = 0.0;
    for (int f : vertex_faces[v]) {
        E += computeFaceEnergy(f) * face_data[f].area;
    }
    return E;
}
// 检查单个面片是否翻转（即雅可比矩阵 J 的行列式是否非正）
bool MIPS::isFaceFlipped(int fid) const {
    Eigen::Matrix2d J = computeJacobian(fid);
    return J.determinant() <= 1e-12;
}

// 数值梯度计算：对面片的每个顶点在 UV 中做微小扰动，计算能量变化来近似梯度
Eigen::Matrix<double, 2, 3> MIPS::numericalFaceGradient(int fid, double eps) const {
    Eigen::Matrix<double, 2, 3> grad;
    grad.setZero();

    int i0 = F(0, fid);
    int i1 = F(1, fid);
    int i2 = F(2, fid);
    int vids[3] = {i0, i1, i2};

    // 因为函数是 const，这里复制一份 UV 来做扰动
    for (int lv = 0; lv < 3; ++lv) {
        int vid = vids[lv];

        for (int d = 0; d < 2; ++d) {
            MIPS temp = *this;

            temp.UV(d, vid) += eps;
            double Ep = temp.computeFaceEnergy(fid);

            temp.UV(d, vid) -= 2.0 * eps;
            double Em = temp.computeFaceEnergy(fid);

            grad(d, lv) = (Ep - Em) / (2.0 * eps);
        }
    }

    return grad;
}

// 计算面的梯度
Eigen::Matrix<double, 2, 3> MIPS::computeFaceGradient(int fid) const {
    return numericalFaceGradient(fid);
}


Eigen::Vector2d MIPS::computeLocalGradient(int v) const {
    Eigen::Vector2d grad(0.0, 0.0);
    for (int f : vertex_faces[v]) {
        Eigen::Matrix<double, 2, 3> fg = computeFaceGradient(f);
        int local_vid = -1;
        for (int lv = 0; lv < 3; ++lv) {
            if (face_data[f].vid[lv] == v) {
                local_vid = lv;
                break;
            }
        }
        if (local_vid != -1) {
            grad += fg.col(local_vid) * face_data[f].area;
        }
    }
    return grad;
}




void MIPS::update_vertex_gd(int v, double init_step) {
    if (is_fixed[v]) return;

    Eigen::Vector2d u_old = UV.col(v);

    Eigen::Vector2d grad = computeLocalGradient(v);
    double gnorm = grad.norm();
    if (gnorm < 1e-12) return;

    double step = init_step;
    double E0 = computeLocalEnergy(v);

    for (int t = 0; t < 20; ++t) {

        UV.col(v) = u_old - step * grad;

        //只检查一环
        bool flipped = false;
        for (int f : vertex_faces[v]) {
            if (isFaceFlipped(f)) {
                flipped = true;
                break;
            }
        }

        if (!flipped) {
            double E1 = computeLocalEnergy(v);
            if (E1 < E0) {
                return; // 成功
            }
        }

        step *= 0.85;   
    }

    UV.col(v) = u_old;
}



bool MIPS::solve(int stages,int inner_iters,
                            double step_size,
                            double tol) {
    // 记录约束初始位置（保险起见）
    for (int i = 0; i < nV; ++i) {
        constraint_initial_uv[i] = UV.col(i);
    }

    constraint_progress = 0.0;
    setConstraintProgress(constraint_progress);
    applyConstraints();

    // ===== 阶段推进 =====
    for (int s = 0; s < stages; ++s) {

        double t = (s + 1.0) / stages;
        advanceConstraintProgressSafely(t, 20, 0.5);

        std::cout << "\n[Stage " << s+1 << "/" << stages << "]\n";

        for (int it = 0; it < inner_iters; ++it) {

            double E0 = computeTotalEnergy();

            for (int c = 0; c < num_colors; ++c) {

                for (int v : color_groups[c]) {
                    update_vertex_gd(v, step_size);
                }
            }

            applyConstraints();

            double E1 = computeTotalEnergy();

            std::cout << "  iter " << it
                    << "  E = " << E1
                    << std::endl;

            if (std::abs(E1 - E0) < tol) break;
        }
    }


    return true;
}