#include "mips.h"

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

// 检查单个面片是否翻转（即雅可比矩阵 J 的行列式是否非正）
bool MIPS::isFaceFlipped(int fid) const {
    Eigen::Matrix2d J = computeJacobian(fid);
    return J.determinant() <= 1e-12;
}

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


Eigen::Matrix2d MIPS::project_to_valid_distortion(const Eigen::Matrix2d& J)
{
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Matrix2d U = svd.matrixU();
    Eigen::Matrix2d V = svd.matrixV();
    Eigen::Vector2d S = svd.singularValues();

    double sigma1 = S(0);
    double sigma2 = S(1);

    // ===== 1. 防止退化 =====
    const double eps = 1e-6;
    sigma1 = std::max(sigma1, eps);
    sigma2 = std::max(sigma2, eps);

    // ===== 2. 防止翻转 =====
    if ((U * V.transpose()).determinant() < 0)
    {
        // 修正 orientation
        V.col(1) *= -1;
        sigma2 *= -1;
    }

    // ===== 3. MIPS-like clamp（关键）=====
    double ratio = sigma1 / sigma2;

    const double max_ratio = 10.0; // 可调

    if (ratio > max_ratio)
    {
        sigma1 = sigma2 * max_ratio;
    }
    if (1.0 / ratio > max_ratio)
    {
        sigma2 = sigma1 * max_ratio;
    }

    // ===== 4. 重建 =====
    Eigen::Matrix2d Sigma;
    Sigma << sigma1, 0,
             0, sigma2;

    return U * Sigma * V.transpose();
}

void MIPS::localStep(){
    local_J.clear();
    local_J.resize(nF);
    for (int f = 0; f < nF; ++f) {
        Eigen::Matrix2d J = computeJacobian(f);
        Eigen::Matrix2d L = project_to_valid_distortion(J);
        local_J[f] = L;
    }
}

bool MIPS::globalStep()
{
    using T = Eigen::Triplet<double>;

    std::vector<T> triplets;
    Eigen::VectorXd b(2 * nV);
    b.setZero();

    // ===== 遍历每个三角形 =====
    for (int f = 0; f < nF; ++f)
    {
        const auto& fd = face_data[f];
        if (!fd.valid) continue;

        int i0 = fd.vid(0);
        int i1 = fd.vid(1);
        int i2 = fd.vid(2);

        double w = fd.area;

        // local step结果
        const Eigen::Matrix2d& L = local_J[f];

        // T = L * Dm
        Eigen::Matrix2d Tmat = L * fd.Dm;

        Eigen::Vector2d t1 = Tmat.col(0);
        Eigen::Vector2d t2 = Tmat.col(1);

        // ===== 约束1: u1 - u0 = t1 =====
        for (int d = 0; d < 2; ++d)
        {
            int row1 = 2 * i1 + d;
            int row0 = 2 * i0 + d;

            triplets.emplace_back(row1, row1,  w);
            triplets.emplace_back(row1, row0, -w);

            triplets.emplace_back(row0, row1, -w);
            triplets.emplace_back(row0, row0,  w);

            b(row1) += w * t1(d);
            b(row0) -= w * t1(d);
        }

        // ===== 约束2: u2 - u0 = t2 =====
        for (int d = 0; d < 2; ++d)
        {
            int row2 = 2 * i2 + d;
            int row0 = 2 * i0 + d;

            triplets.emplace_back(row2, row2,  w);
            triplets.emplace_back(row2, row0, -w);

            triplets.emplace_back(row0, row2, -w);
            triplets.emplace_back(row0, row0,  w);

            b(row2) += w * t2(d);
            b(row0) -= w * t2(d);
        }
    }

    // ===== 构建稀疏矩阵 =====
    Eigen::SparseMatrix<double> A(2 * nV, 2 * nV);
    A.setFromTriplets(triplets.begin(), triplets.end());

    // ===== 加入硬约束（非常关键）=====
    const double w_fix = 1e8;

    for (int i = 0; i < nV; ++i)
    {
        if (is_fixed[i])
        {
            for (int d = 0; d < 2; ++d)
            {
                int idx = 2 * i + d;

                A.coeffRef(idx, idx) += w_fix;
                b(idx) += w_fix * constraint_active_uv[i](d);
                
            }
        }
    }

    // ===== 求解 =====
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);

    if (solver.info() != Eigen::Success)
    {
        std::cout << "Decomposition failed\n";
        return false;
    }

    Eigen::VectorXd x = solver.solve(b);

    if (solver.info() != Eigen::Success)
    {
        std::cout << "Solve failed\n";
        return false;
    }

    // ===== 写回 UV =====
    for (int i = 0; i < nV; ++i)
    {
        UV(0, i) = x(2 * i);
        UV(1, i) = x(2 * i + 1);
    }

    return true;
}

Eigen::Matrix2Xd MIPS::AndersonAcceleration(
    const Eigen::Matrix2Xd& UV_k,
    const Eigen::Matrix2Xd& G_k,
    const Eigen::Matrix2Xd& F_k)
{
    (void)UV_k;
    int k = static_cast<int>(F_history.size());

    // ===== 如果历史不够，直接返回 LG =====
    if (k < 2)
        return G_k;

    // k 条历史只能构造 k-1 条差分，避免访问到 [-1]
    int m = std::min(k - 1, anderson_history_m);
    if (m <= 0)
        return G_k;

    int dim = 2 * nV;

    // ===== 构造 ΔF =====
    Eigen::MatrixXd D(dim, m);

    for (int j = 0; j < m; ++j)
    {
        const auto& F1 = F_history[k - j - 1];
        const auto& F0 = F_history[k - j - 2];

        if (F1.rows() != F0.rows() || F1.cols() != F0.cols()) {
            return G_k;
        }

        Eigen::Matrix2Xd dF_mat = F1 - F0;
        D.col(j) = Eigen::Map<const Eigen::VectorXd>(dF_mat.data(), dim);
    }

    // 当前 residual
    Eigen::VectorXd fk = Eigen::Map<const Eigen::VectorXd>(
        F_k.data(), dim);

    // ===== 解最小二乘：min || fk - D θ ||^2 =====
    Eigen::VectorXd theta;

    // 正规方程（m很小，没问题）
    Eigen::MatrixXd DtD = D.transpose() * D;
    Eigen::VectorXd Dtf = D.transpose() * fk;

    // 防止奇异（很关键）
    double reg = 1e-8;
    DtD += reg * Eigen::MatrixXd::Identity(m, m);

    theta = DtD.ldlt().solve(Dtf);

    // ===== 构造 ΔG =====
    Eigen::MatrixXd DG(dim, m);

    for (int j = 0; j < m; ++j)
    {
        const auto& G1 = G_history[k - j - 1];
        const auto& G0 = G_history[k - j - 2];

        if (G1.rows() != G0.rows() || G1.cols() != G0.cols()) {
            return G_k;
        }

        Eigen::Matrix2Xd dG_mat = G1 - G0;
        DG.col(j) = Eigen::Map<const Eigen::VectorXd>(dG_mat.data(), dim);
    }

    // ===== Anderson 更新 =====
    Eigen::VectorXd gk = Eigen::Map<const Eigen::VectorXd>(
        G_k.data(), dim);

    Eigen::VectorXd x_new = gk - DG * theta;

    // reshape 回 UV
    Eigen::Matrix2Xd UV_new(2, nV);
    UV_new = Eigen::Map<Eigen::Matrix<double, 2, Eigen::Dynamic>>(x_new.data(), 2, nV);

    return UV_new;
}

bool MIPS::solve(int stages,int inner_iters,
                            int final_iters,
                            double step_size,
                            double tol) {
    // 记录约束初始位置（保险起见）
    for (int i = 0; i < nV; ++i) {
        constraint_initial_uv[i] = UV.col(i);
    }
    G_history.clear();
    F_history.clear();
    Q_history.clear();

    constraint_progress = 0.0;
    setConstraintProgress(constraint_progress);
    applyConstraints();

    for (int s = 0; s < stages; ++s)
    {
        double target_t = double(s + 1) / stages;
        advanceConstraintProgressSafely(target_t, 20, 0.5);

        std::cout << "\n[Stage " << s + 1 << "/" << stages << "]\n";

        for (int it = 0; it < inner_iters; ++it)
        {
            Eigen::Matrix2Xd UV_old = UV;

            // ===== local-global =====
            localStep();
            globalStep();
            Eigen::Matrix2Xd UV_LG = UV;

            // ===== residual =====
            Eigen::Matrix2Xd F = UV_LG - UV_old;
            double res_norm = F.norm();

            double E_lg = computeTotalEnergy();

            std::cout << "  iter " << it
                    << "  E = " << E_lg
                    << "  |F| = " << res_norm << std::endl;

            if (res_norm < tol)
                break;

            // ===== Anderson =====
            Eigen::Matrix2Xd UV_new = UV_LG;

            if (F_history.size() >= 2)   // ⚠️ 必须 ≥2
            {
                Eigen::Matrix2Xd UV_AA = AndersonAcceleration(UV_old, UV_LG, F);

                // ===== energy check =====
                UV = UV_AA;
                applyConstraints();

                double E_aa = computeTotalEnergy();

                if (E_aa < E_lg && !hasAnyFlippedFace())
                {
                    UV_new = UV;   // 接受 AA
                }
                else
                {
                    UV_new = UV_LG; // 回退 LG
                }
            }

            // ===== 写回最终解 =====
            UV = UV_new;
            applyConstraints();

            // ===== 存历史 =====
            Q_history.push_back(UV_old);
            G_history.push_back(UV_LG);
            F_history.push_back(F);

            if (F_history.size() > anderson_history_m + 1)
            {
                Q_history.erase(Q_history.begin());
                G_history.erase(G_history.begin());
                F_history.erase(F_history.begin());
            }
        }
    }

    // ===== 最终 refinement =====
    std::cout << "\n[Final refinement]\n";

    double prevE = computeTotalEnergy();

    for (int it = 0; it < final_iters; ++it)
    {
        Eigen::Matrix2Xd UV_old = UV;

        localStep();
        globalStep();

        Eigen::Matrix2Xd UV_LG = UV;
        Eigen::Matrix2Xd F = UV_LG - UV_old;

        double res = F.norm();
        double E_lg = computeTotalEnergy();

        std::cout << "  final " << it
                  << "  E = " << E_lg
                  << "  |F| = " << res << std::endl;

        if (res < tol)
        {
            std::cout << "Converged by residual.\n";
            return true;
        }

        // ===== Anderson =====
        Eigen::Matrix2Xd UV_AA = UV_LG;

        if (F_history.size() >= 1)
        {
            UV_AA = AndersonAcceleration(UV_old, UV_LG, F);
        }

        // energy check
        UV = UV_AA;
        applyConstraints();

        double E_aa = computeTotalEnergy();

        if (E_aa < E_lg && !hasAnyFlippedFace())
        {
            // accept
        }
        else
        {
            UV = UV_LG;
            applyConstraints();
        }

        // 收敛判断
        double E = computeTotalEnergy();
        if (std::abs(prevE - E) < tol)
        {
            std::cout << "Converged by energy change.\n";
            return true;
        }

        prevE = E;

        // 存历史
        Q_history.push_back(UV_old);
        G_history.push_back(UV_LG);
        F_history.push_back(F);

        if (F_history.size() > anderson_history_m + 1)
        {
            Q_history.erase(Q_history.begin());
            G_history.erase(G_history.begin());
            F_history.erase(F_history.begin());
        }
    }

    return true;
}