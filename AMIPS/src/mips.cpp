#include "mips.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

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
    buildOneRingAdjacency();
    computeVertexColors();
    total_energy = computeTotalEnergy();
}


bool MIPS::buildLocal2DReference(int dim) {
    face_data.resize(nF);
    if(dim != 2 && dim != 3) {
        std::cerr << "Unsupported dimension for local reference: " << dim << std::endl;
        return false;
    }
    if (dim == 2) {
        for (int f = 0; f < nF; ++f) {
            int i = F(0, f);
            int j = F(1, f);
            int k = F(2, f);

            Eigen::Vector2d p0, p1, p2;
            p0 = V.col(i).head<2>();
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

void MIPS::buildOneRingAdjacency() {
    vertex_to_faces.assign(nV, std::vector<int>());
    for (int f = 0; f < nF; ++f) {
        vertex_to_faces[F(0, f)].push_back(f);
        vertex_to_faces[F(1, f)].push_back(f);
        vertex_to_faces[F(2, f)].push_back(f);
    }
}

void MIPS::computeVertexColors() {
    vertex_colors.assign(nV, -1);
    const int MAX_COLORS = 4;

    for (int vid = 0; vid < nV; ++vid) {
        std::vector<bool> used_color(MAX_COLORS, false);
        const auto& faces = vertex_to_faces[vid];
        for (int fid : faces) {
            const Eigen::Vector3i& tri = face_data[fid].vid;
            for (int k = 0; k < 3; ++k) {
                int neighbor = tri[k];
                if (neighbor != vid && vertex_colors[neighbor] >= 0) {
                    if (vertex_colors[neighbor] < MAX_COLORS) {
                        used_color[vertex_colors[neighbor]] = true;
                    }
                }
            }
        }

        for (int c = 0; c < MAX_COLORS; ++c) {
            if (!used_color[c]) {
                vertex_colors[vid] = c;
                break;
            }
        }

        if (vertex_colors[vid] < 0) {
            vertex_colors[vid] = vid % MAX_COLORS;
        }
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

double MIPS::computeFaceMIPS(int fid) const {
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

double MIPS::computeFaceAMIPS(int fid) const {
    if (!face_data[fid].valid) {
        return 0.0;
    }

    Eigen::Matrix2d J = computeJacobian(fid);
    double detJ = J.determinant();
    if (detJ <= 1e-12) {
        return 1e20;
    }

    Eigen::JacobiSVD<Eigen::Matrix2d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector2d sigma = svd.singularValues();
    double sigma1 = sigma(0);
    double sigma2 = sigma(1);

    if (sigma2 < 1e-12) {
        return 1e20;
    }

    return (sigma1 / sigma2) + (sigma2 / sigma1);
}

double computeDetEnergy(double detJ) {
    if (detJ <= 1e-12) {
        return 1e20;
    }
    return 0.5 * (detJ + 1.0 / detJ);
}

double MIPS::computeFaceIsometric(int fid) const {
    if (!face_data[fid].valid) {
        return 0.0;
    }

    Eigen::Matrix2d J = computeJacobian(fid);
    double detJ = J.determinant();
    if (detJ <= 1e-12) {
        return 1e20;
    }

    double E_mips = computeFaceMIPS(fid);
    double E_det = computeDetEnergy(detJ);
    return iso_alpha * E_mips + (1.0 - iso_alpha) * E_det;
}

double MIPS::computeFaceEnergy(int fid) const {
    if (!face_data[fid].valid) {
        return 0.0;
    }

    if (energy_type == EnergyType::MIPS) {
        return computeFaceMIPS(fid);
    } else if (energy_type == EnergyType::AMIPS) {
        return computeFaceAMIPS(fid);
    } else {
        return computeFaceIsometric(fid);
    }
}

double MIPS::computeAMIPSWeight(int fid) const {
    if (!face_data[fid].valid || penalty_scale <= 0.0) {
        return 1.0;
    }

    double Ef = computeFaceEnergy(fid);
    if (Ef > 1e19) {
        return 1e20;
    }

    return std::exp(penalty_scale * Ef);
}

static constexpr double BARRIER_EPS = 1e-6;

static double computeFaceLogBarrier(double detJ, double lambda) {
    if (lambda <= 0.0) return 0.0;
    if (detJ <= BARRIER_EPS) {
        return lambda * (-std::log(BARRIER_EPS));
    }
    return lambda * (-std::log(detJ + BARRIER_EPS));
}

double MIPS::computeTotalEnergy() const {
    double total = 0.0;

    for (int f = 0; f < nF; ++f) {
        double Ef = computeFaceEnergy(f);
        if (Ef > 1e10) continue;

        double weight = 1.0;
        if (penalty_scale > 0.0 && energy_type == EnergyType::AMIPS) {
            weight = computeAMIPSWeight(f);
        }

        total += face_data[f].area * weight * Ef;

        if (log_barrier_lambda > 0.0 && face_data[f].valid) {
            Eigen::Matrix2d J = computeJacobian(f);
            total += face_data[f].area * computeFaceLogBarrier(J.determinant(), log_barrier_lambda);
        }
    }

    return total;
}

bool MIPS::isFaceFlipped(int fid) const {
    Eigen::Matrix2d J = computeJacobian(fid);
    return J.determinant() <= 1e-12;
}

// =============================================================================
// GRADIENT COMPUTATION - CRITICAL FIX
//
// Paper defines: E* = Σ exp(s * E_f)
//
// CORRECT chain rule:
//   ∂(exp(sE)) / ∂x = exp(sE) * ∂(sE) / ∂x = s * exp(sE) * ∂E / ∂x
//
// The gradient contribution from face f to vertex v is:
//   w_f = exp(s * E_f)           (AMIPS weight)
//   ∇_v E_f = base gradient
//   ∇_v (w_f * E_f) = w_f * s * ∇_v E_f   (chain rule!)
//
// Fix: multiply by s (penalty_scale) in addition to w_f
// =============================================================================

Eigen::Matrix<double, 2, 3> MIPS::computeFaceBaseGradient(int fid) {
    Eigen::Matrix<double, 2, 3> grad;
    grad.setZero();

    const auto& tri = face_data[fid].vid;
    int i0 = tri[0], i1 = tri[1], i2 = tri[2];
    int vids[3] = {i0, i1, i2};

    constexpr double eps = 1e-6;

    for (int lv = 0; lv < 3; ++lv) {
        int vid = vids[lv];
        for (int d = 0; d < 2; ++d) {
            double orig = UV(d, vid);
            UV(d, vid) = orig + eps;
            double Ep = computeFaceEnergy(fid);

            UV(d, vid) = orig - eps;
            double Em = computeFaceEnergy(fid);

            UV(d, vid) = orig;

            grad(d, lv) = (Ep - Em) / (2.0 * eps);
        }
    }

    return grad;
}

Eigen::Matrix<double, 2, 3> MIPS::computeFaceGradient(int fid) {
    if (!face_data[fid].valid) {
        return Eigen::Matrix<double, 2, 3>::Zero();
    }

    Eigen::Matrix<double, 2, 3> grad = computeFaceBaseGradient(fid);

    double weight = 1.0;
    double s = penalty_scale;

    if (penalty_scale > 0.0 && energy_type == EnergyType::AMIPS) {
        weight = computeAMIPSWeight(fid);

        // FIX: multiply by s for chain rule
        // ∂(exp(sE))/∂x = s * exp(sE) * ∂E/∂x = s * weight * ∂E/∂x
        grad = s * weight * grad;
    } else {
        grad *= face_data[fid].area;
    }

    // Optional log-barrier gradient (NOT in paper)
    if (log_barrier_lambda > 0.0 && face_data[fid].valid) {
        constexpr double eps = 1e-6;
        const Eigen::Vector3i& tri = face_data[fid].vid;

        Eigen::Matrix2d J = computeJacobian(fid);
        double detJ = J.determinant();
        double denom = detJ + BARRIER_EPS;
        double barrier_weight = log_barrier_lambda * face_data[fid].area / denom;

        for (int lv = 0; lv < 3; ++lv) {
            int v_local = tri[lv];

            for (int dim = 0; dim < 2; ++dim) {
                double orig = UV(dim, v_local);
                UV(dim, v_local) = orig + eps;
                double detJp = computeJacobian(fid).determinant();

                UV(dim, v_local) = orig - eps;
                double detJm = computeJacobian(fid).determinant();

                UV(dim, v_local) = orig;

                double d_determinant = (detJp - detJm) / (2.0 * eps);
                grad(dim, lv) += -barrier_weight * d_determinant;
            }
        }
    }

    return grad;
}

Eigen::Matrix<double, 2, 3> MIPS::numericalFaceGradient(int fid, double eps) const {
    Eigen::Matrix<double, 2, 3> grad;
    grad.setZero();

    int i0 = F(0, fid);
    int i1 = F(1, fid);
    int i2 = F(2, fid);
    int vids[3] = {i0, i1, i2};

    for (int lv = 0; lv < 3; ++lv) {
        int vid = vids[lv];

        for (int d = 0; d < 2; ++d) {
            double orig = UV(d, vid);

            UV(d, vid) = orig + eps;
            double Ep = computeFaceEnergy(fid);

            UV(d, vid) = orig - eps;
            double Em = computeFaceEnergy(fid);

            UV(d, vid) = orig;

            grad(d, lv) = (Ep - Em) / (2.0 * eps);
        }
    }

    return grad;
}

void MIPS::computeGradient(Eigen::Matrix2Xd& grad) const {
    grad.resize(2, nV);
    grad.setZero();

    for (int f = 0; f < nF; ++f) {
        Eigen::Matrix<double, 2, 3> fg = numericalFaceGradient(f);

        int i0 = F(0, f);
        int i1 = F(1, f);
        int i2 = F(2, f);

        grad.col(i0) += fg.col(0);
        grad.col(i1) += fg.col(1);
        grad.col(i2) += fg.col(2);
    }

    for (int i = 0; i < nV; ++i) {
        if (is_fixed[i]) {
            grad.col(i).setZero();
        }
    }
}

Eigen::Vector2d MIPS::computeVertexGradient(int vid) const {
    Eigen::Vector2d g = Eigen::Vector2d::Zero();

    const auto& faces = vertex_to_faces[vid];
    for (int fid : faces) {
        Eigen::Matrix<double, 2, 3> fg = numericalFaceGradient(fid);
        int lv = -1;
        if (F(0, fid) == vid) lv = 0;
        else if (F(1, fid) == vid) lv = 1;
        else if (F(2, fid) == vid) lv = 2;
        if (lv >= 0) {
            g += fg.col(lv);
        }
    }

    return g;
}

// =============================================================================
// Line Search - Paper Sec 3.2
//
// 论文描述:
//   "The initial λ is chosen such that p - λ∇F is on the boundary of the 
//    one ring region of p. λ is reduced by 85% if F increases or one of the
//    neighboring simplices of p inverts."
// =============================================================================

double MIPS::computeVertexOneRingBoundary(int vid) const {
    // 计算顶点 one-ring 的边界长度
    // 论文: "λ 使 p - λ∇F 在 one-ring 边界上"
    double boundary = 0.0;
    const auto& faces = vertex_to_faces[vid];

    for (int fid : faces) {
        const Eigen::Vector3i& tri = face_data[fid].vid;
        for (int k = 0; k < 3; ++k) {
            int vj = tri[k];
            if (vj != vid) {
                boundary += (UV.col(vid) - UV.col(vj)).norm();
            }
        }
    }

    return boundary / static_cast<double>(faces.size());
}

bool MIPS::wouldCauseFlip(int vid, const Eigen::Vector2d& new_pos) const {
    // 检查 one-ring 是否有 flip
    const auto& faces = vertex_to_faces[vid];
    Eigen::Vector2d old_pos = UV.col(vid);

    UV.col(vid) = new_pos;

    for (int fid : faces) {
        if (isFaceFlipped(fid)) {
            UV.col(vid) = old_pos;
            return true;
        }
    }

    UV.col(vid) = old_pos;
    return false;
}

// 论文 Sec 3.2: 逐顶点 line search
bool MIPS::lineSearchForVertex(int vid, const Eigen::Vector2d& grad, double& step_size) {
    if (is_fixed[vid] || grad.norm() < 1e-12) {
        return false;
    }

    // 论文: "λ such that p - λ∇F is on the boundary of the one ring region"
    double boundary_len = computeVertexOneRingBoundary(vid);
    double lambda = boundary_len;  // 初始 λ = one-ring 边界长度

    const int max_ls_iters = 20;

    double E_before = computeTotalEnergy();
    Eigen::Vector2d old_pos = UV.col(vid);

    for (int ls_iter = 0; ls_iter < max_ls_iters; ++ls_iter) {
        // 论文: p_new = p - λ * ∇F
        Eigen::Vector2d new_pos = old_pos - lambda * grad;

        // 检查 flip (论文: "one of the neighboring simplices of p inverts")
        if (wouldCauseFlip(vid, new_pos)) {
            lambda *= 0.85;  // 论文: "reduced by 85%"
            continue;
        }

        // 应用更新
        UV.col(vid) = new_pos;
        double E_after = computeTotalEnergy();

        // 论文: "if F increases"
        if (E_after <= E_before) {
            step_size = lambda;
            return true;
        }

        // 论文: "λ is reduced by 85% if F increases"
        UV.col(vid) = old_pos;
        lambda *= 0.85;
    }

    return false;
}

// =============================================================================
// BCD MODE 0: COLORED BLOCKS (Paper Sec 3.2 Default Implementation)
//
// 论文:
//   "We partition the mesh vertices into several blocks where any two vertices
//    in the same block have no connected edges in the mesh."
//   "For each vertex p in a block, we update it by one step of gradient descent:
//    p_new := p - λ ∇_p F"
//
// 实现:
//   for each color block:
//       for each vertex in block:
//           compute gradient
//           line search
//       apply all updates for this block simultaneously
// =============================================================================
int MIPS::gradientDescentSweep_GaussSeidel(double /* step_base */) {
    int updated = 0;
    int max_color = 0;
    for (int c : vertex_colors) {
        max_color = std::max(max_color, c);
    }

    // 论文: "For each block" - 按颜色分块处理
    for (int c = 0; c <= max_color; ++c) {
        // 收集当前颜色的所有顶点
        std::vector<int> block_vertices;
        for (int vid = 0; vid < nV; ++vid) {
            if (vertex_colors[vid] == c && !is_fixed[vid]) {
                block_vertices.push_back(vid);
            }
        }

        if (block_vertices.empty()) continue;

        // 论文: "For each vertex p in a block"
        // 计算每个顶点的梯度
        std::vector<Eigen::Vector2d> grads(block_vertices.size());
        for (size_t i = 0; i < block_vertices.size(); ++i) {
            grads[i] = computeVertexGradient(block_vertices[i]);
        }

        // 逐顶点 line search
        std::vector<Eigen::Vector2d> deltas(block_vertices.size());
        for (size_t i = 0; i < block_vertices.size(); ++i) {
            int vid = block_vertices[i];
            double lambda_used;
            if (lineSearchForVertex(vid, grads[i], lambda_used)) {
                deltas[i] = -lambda_used * grads[i];
                ++updated;
            } else {
                deltas[i].setZero();
            }
        }

        // 论文: 同一 block 内顶点可同时更新（不相邻所以独立）
        for (size_t i = 0; i < block_vertices.size(); ++i) {
            if (deltas[i].norm() > 1e-14) {
                UV.col(block_vertices[i]) += deltas[i];
            }
        }
    }

    applyConstraints();
    return updated;
}

// 保留 Jacobi 版本作为可选快速版本（论文不推荐但有时更快）
int MIPS::gradientDescentSweep_Jacobi(double step_base) {
    std::vector<Eigen::Vector2d> grad_snap(nV);
    for (int vid = 0; vid < nV; ++vid) {
        if (is_fixed[vid]) {
            grad_snap[vid].setZero();
        } else {
            grad_snap[vid] = computeVertexGradient(vid);
        }
    }

    std::vector<Eigen::Vector2d> delta(nV);
    int movable = 0;

    for (int vid = 0; vid < nV; ++vid) {
        if (is_fixed[vid]) {
            delta[vid].setZero();
            continue;
        }

        const Eigen::Vector2d& g = grad_snap[vid];
        if (g.norm() < 1e-12) {
            delta[vid].setZero();
            continue;
        }

        double boundary_len = computeVertexOneRingBoundary(vid);
        delta[vid] = -step_base * boundary_len * g;
        ++movable;
    }

    if (movable == 0) return 0;

    for (int vid = 0; vid < nV; ++vid) {
        if (delta[vid].norm() > 1e-14) {
            UV.col(vid) += delta[vid];
        }
    }
    applyConstraints();

    double E_new = computeTotalEnergy();

    if (E_new > E_last_sweep * 50.0) {
        for (int vid = 0; vid < nV; ++vid) {
            if (delta[vid].norm() > 1e-14) {
                UV.col(vid) -= delta[vid];
            }
        }
        applyConstraints();
        return 0;
    }

    E_last_sweep = E_new;
    return movable;
}

// =============================================================================
// BCD MODE 2: COLORED BLOCKS (Parallel-friendly)
// =============================================================================
int MIPS::gradientDescentSweep_Colored(double step_base) {
    int updated = 0;
    int max_color = 0;
    for (int c : vertex_colors) {
        max_color = std::max(max_color, c);
    }

    for (int c = 0; c <= max_color; ++c) {
        std::vector<int> color_vertices;
        for (int vid = 0; vid < nV; ++vid) {
            if (vertex_colors[vid] == c && !is_fixed[vid]) {
                color_vertices.push_back(vid);
            }
        }

        if (color_vertices.empty()) continue;

        std::vector<Eigen::Vector2d> grads(color_vertices.size());
        for (size_t i = 0; i < color_vertices.size(); ++i) {
            grads[i] = computeVertexGradient(color_vertices[i]);
        }

        std::vector<Eigen::Vector2d> updates(color_vertices.size());
        for (size_t i = 0; i < color_vertices.size(); ++i) {
            int vid = color_vertices[i];
            double step_size = step_base;
            if (lineSearchForVertex(vid, grads[i], step_size)) {
                updates[i] = -step_size * grads[i];
                ++updated;
            } else {
                updates[i].setZero();
            }
        }

        for (size_t i = 0; i < color_vertices.size(); ++i) {
            if (updates[i].norm() > 1e-14) {
                UV.col(color_vertices[i]) += updates[i];
            }
        }
    }

    applyConstraints();
    return updated;
}

int MIPS::gradientDescentSweep(double step_base, int mode) {
    switch (mode) {
        case 0: return gradientDescentSweep_GaussSeidel(step_base);
        case 1: return gradientDescentSweep_Jacobi(step_base);
        case 2: return gradientDescentSweep_Colored(step_base);
        default: return gradientDescentSweep_GaussSeidel(step_base);
    }
}

int MIPS::gradientDescentSweep(double step_base) {
    return gradientDescentSweep_GaussSeidel(step_base);
}

bool MIPS::advanceConstraintProgressSoft(double target_t,
                                       double E_bound,
                                       int repair_iters,
                                       int max_attempts) {
    target_t = std::max(constraint_progress, std::min(1.0, target_t));
    if (target_t <= constraint_progress) return true;

    double cand_t = target_t;
    double shrink = 0.5;

    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        setConstraintProgress(cand_t);
        applyConstraints();

        double E_before = computeTotalEnergy();
        for (int r = 0; r < repair_iters; ++r) {
            int updated = gradientDescentSweep(1e-2);
            if (updated == 0) break;
        }
        double E_after = computeTotalEnergy();

        if (E_after < E_bound) {
            constraint_progress = cand_t;
            std::cout << "  soft advance: t=" << std::fixed << std::setprecision(4) << cand_t
                      << " E_before=" << E_before << " E_after=" << E_after
                      << " (bound=" << E_bound << ")" << std::endl;
            return true;
        }

        setConstraintProgress(constraint_progress);
        applyConstraints();

        double remaining = target_t - constraint_progress;
        if (remaining < 1e-8) break;

        cand_t = constraint_progress + remaining * shrink;
        std::cout << "  soft advance: t=" << std::fixed << std::setprecision(4) << cand_t
                  << " energy exploded (" << E_after << " > " << E_bound
                  << "), retrying with smaller step" << std::endl;
    }

    setConstraintProgress(constraint_progress);
    applyConstraints();
    return false;
}

bool MIPS::solve(int stages, int inner_iters,
                 int final_iters,
                 double step_size,
                 double tol,
                 int update_mode) {
    for (int i = 0; i < nV; ++i) {
        constraint_initial_uv[i] = UV.col(i);
    }

    constraint_progress = 0.0;
    setConstraintProgress(constraint_progress);
    applyConstraints();

    double E_ref = computeTotalEnergy();
    double E_bound = std::max(E_ref * 50.0, 1e8);
    E_last_sweep = E_ref;

    std::cout << "=== MIPS Solver ===" << std::endl;
    std::cout << "AMIPS variant: " << (amips_variant == AMIPSVariant::STANDARD ? "STANDARD (paper)" : "BARRIER (comparison)") << std::endl;
    std::cout << "Update mode: " << update_mode << " ";
    if (update_mode == 0) std::cout << "(Gauss-Seidel BCD)";
    else if (update_mode == 1) std::cout << "(Jacobi BCD)";
    else if (update_mode == 2) std::cout << "(Colored blocks)";
    std::cout << std::endl;
    std::cout << "penalty_scale (s): " << penalty_scale << std::endl;
    std::cout << "log_barrier_lambda: " << log_barrier_lambda << std::endl;
    std::cout << "Line search: shrink=" << ls_shrink_factor << std::endl;
    std::cout << "Initial energy: " << E_ref << std::endl;

    for (int s = 0; s < stages; ++s) {
        double t = static_cast<double>(s + 1) / static_cast<double>(stages);

        bool stage_ok = advanceConstraintProgressSoft(t, E_bound, 20, 10);
        if (!stage_ok) {
            std::cout << "  warning: cannot advance to t=" << t
                      << ", staying at t=" << constraint_progress << std::endl;
        }

        std::cout << "\n[Stage " << (s + 1) << "/" << stages
                  << "  t=" << constraint_progress << "]\n";

        double E0 = computeTotalEnergy();
        E_last_sweep = E0;

        for (int it = 0; it < inner_iters; ++it) {
            double E_before = computeTotalEnergy();
            std::cout << "  inner " << it << "  E = " << E_before << std::endl;

            int updated = gradientDescentSweep(step_size);

            double E_after = computeTotalEnergy();

            if (updated == 0) {
                double retry = gradientDescentSweep(step_size * 0.5);
                if (retry == 0) {
                    std::cout << "  no progress, stopping." << std::endl;
                    break;
                }
                E_after = computeTotalEnergy();
            }

            double rel_change = std::abs(E_after - E_before)
                              / std::max(std::abs(E_before), 1.0);
            if (rel_change < tol) {
                std::cout << "  Converged rel=" << rel_change << std::endl;
                break;
            }

            if (E_after > E_bound) {
                std::cout << "  Energy exploded (" << E_after
                          << " > " << E_bound << "), stopping stage." << std::endl;
                break;
            }

            if (E_after > E_before * 1.01) {
                step_size *= 0.8;
                std::cout << "  Energy rose slightly, reducing step to " << step_size << std::endl;
            }
        }

        double E_stage_end = computeTotalEnergy();
        E_bound = std::max(E_stage_end * 10.0, E_bound);
        E_last_sweep = E_stage_end;
    }

    if (!advanceConstraintProgressSoft(1.0, E_bound, 30, 15)) {
        std::cout << "warning: cannot reach full target, "
                  << "final progress = " << constraint_progress << std::endl;
    }

    std::cout << "\n[Final refinement]\n";

    double prevE = computeTotalEnergy();
    E_last_sweep = prevE;
    for (int it = 0; it < final_iters; ++it) {
        double E0 = computeTotalEnergy();
        std::cout << "  final " << it << "  E = " << E0 << std::endl;

        int updated = gradientDescentSweep(step_size);

        double E1 = computeTotalEnergy();
        double rel_change = std::abs(E1 - E0) / std::max(std::abs(E0), 1.0);
        if (rel_change < tol) {
            std::cout << "Final converged rel=" << rel_change << std::endl;
            break;
        }
        if (updated == 0) {
            double retry = gradientDescentSweep(step_size * 0.5);
            if (retry == 0) break;
        }

        prevE = E1;
    }

    total_energy = computeTotalEnergy();
    std::cout << "Done. Final energy = " << total_energy << std::endl;
    std::cout << "Flipped faces: " << (hasAnyFlippedFace() ? "YES (bad!)" : "NO (good!)") << std::endl;
    return true;
}
