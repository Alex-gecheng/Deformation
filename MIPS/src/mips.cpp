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

    buildLocal2DReference(3);
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

void MIPS::applyConstraints() {
    for (const auto& c : constraints) {
        int vid = c.first;
        if (vid >= 0 && vid < nV) {
            UV.col(vid) = constraint_active_uv[vid];
        }
    }
}

bool MIPS::hasAnyFlippedFace() const {
    for (int f = 0; f < nF; ++f) {
        if (isFaceFlipped(f)) {
            return true;
        }
    }
    return false;
}

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

double MIPS::computeTotalEnergy() const {
    double total = 0.0;
    for (int f = 0; f < nF; ++f) {
        total += computeFaceEnergy(f)*face_data[f].area;    //乘以面积权重
    }
    return total;
}

bool MIPS::isFaceFlipped(int fid) const {
    Eigen::Matrix2d J = computeJacobian(fid);
    return J.determinant() <= 1e-12;
}

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

Eigen::Matrix<double, 2, 3> MIPS::computeFaceGradient(int fid) const {
    return numericalFaceGradient(fid);
}

void MIPS::computeGradient(Eigen::Matrix2Xd& grad) const {
    grad.resize(2, nV);
    grad.setZero();

    for (int f = 0; f < nF; ++f) {
        Eigen::Matrix<double, 2, 3> fg = computeFaceGradient(f);
        fg *= face_data[f].area;

        int i0 = F(0, f);
        int i1 = F(1, f);
        int i2 = F(2, f);

        grad.col(i0) += fg.col(0);
        grad.col(i1) += fg.col(1);
        grad.col(i2) += fg.col(2);
    }

    // 约束点梯度清零（不允许更新）
    for (int i = 0; i < nV; ++i) {
        if (is_fixed[i]) {
            grad.col(i).setZero();
        }
    }
}

double MIPS::backtrackingLineSearch(const Eigen::Matrix2Xd& grad,
                                    double init_step,
                                    double shrink,
                                    int max_trials) {
    double E0 = computeTotalEnergy();
    Eigen::Matrix2Xd UV_old = UV;

    double step = init_step;

    for (int t = 0; t < max_trials; ++t) {
        UV = UV_old - step * grad;
        applyConstraints();

        bool flipped = false;
        for (int f = 0; f < nF; ++f) {
            if (isFaceFlipped(f)) {
                flipped = true;
                break;
            }
        }

        if (!flipped) {
            double E1 = computeTotalEnergy();
            if (E1 < E0) {
                return step;
            }
        }

        step *= shrink;
    }

    UV = UV_old;
    return 0.0;
}

bool MIPS::gradientDescentStep(double step_size) {
    Eigen::Matrix2Xd grad;
    computeGradient(grad);

    double gnorm = grad.norm();
    if (gnorm < 1e-12) return false;

    Eigen::Matrix2Xd UV_old = UV;
    double step = backtrackingLineSearch(grad, step_size);

    if (step <= 0.0) {
        UV = UV_old;
        return false;
    }

    return true;
}

bool MIPS::solve(int stages,int inner_iters,
                            int final_iters,
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
        // 尝试安全推进到当前阶段目标，若翻转则自动回退步长。
        double t = static_cast<double>(s + 1) / static_cast<double>(stages);
        bool stage_ok = advanceConstraintProgressSafely(t, 20, 0.5);
        if (!stage_ok) {
            std::cout << "  warning: cannot safely reach stage target t=" << t
                      << ", stay at t=" << constraint_progress << std::endl;
        }

        std::cout << "\n[Stage " << (s + 1) << "/" << stages << "]\n";

        // 每阶段做多次优化
        for (int it = 0; it < inner_iters; ++it) {
            Eigen::Matrix2Xd grad;
            computeGradient(grad);

            double gnorm = grad.norm();
            double E0 = computeTotalEnergy();

            std::cout << "  inner " << it
                      << "  E = " << E0
                      << "  |grad| = " << gnorm << std::endl;

            if (gnorm < tol) break;

            double step = backtrackingLineSearch(grad, step_size, 0.5, 30);
            if (step <= 0.0) {
                std::cout << "  line search failed at stage " << s + 1 << "\n";
                return false;
            }

            // 注意：你当前 backtrackingLineSearch 已经更新了 UV
            applyConstraints();
        }
    }

    // ===== 最终安全推进到目标后，完整优化 =====
    if (!advanceConstraintProgressSafely(1.0, 40, 0.5)) {
        std::cout << "warning: cannot safely reach full target constraints, "
                  << "final progress = " << constraint_progress << std::endl;
    }

    std::cout << "\n[Final refinement]\n";

    double prevE = computeTotalEnergy();
    for (int it = 0; it < final_iters; ++it) {
        Eigen::Matrix2Xd grad;
        computeGradient(grad);

        double gnorm = grad.norm();
        std::cout << "  final " << it
                  << "  E = " << prevE
                  << "  |grad| = " << gnorm << std::endl;

        if (gnorm < tol) {
            std::cout << "Converged by gradient norm.\n";
            return true;
        }

        double step = backtrackingLineSearch(grad, step_size, 0.5, 30);
        if (step <= 0.0) {
            std::cout << "Final line search failed.\n";
            return false;
        }

        double E = computeTotalEnergy();
        if (std::abs(prevE - E) < tol) {
            std::cout << "Converged by energy change.\n";
            return true;
        }

        prevE = E;
        applyConstraints();
    }

    return true;
}