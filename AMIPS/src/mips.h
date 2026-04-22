#ifndef MIPS_H
#define MIPS_H

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <limits>

class MIPS {
public:
    using Constraint = std::pair<int, Eigen::Vector3d>;

    enum class EnergyType { MIPS, AMIPS, ISOMETRIC };

    enum class AMIPSVariant { STANDARD, BARRIER };

    MIPS(const Eigen::Matrix3Xd& V,
        const Eigen::Matrix3Xi& F,
        const std::vector<Constraint>& constraints);

    ~MIPS() = default;

    bool solve(int stages, int inner_iters, int final_iters,
               double step_size = 1e-2,
               double tol = 1e-8,
               int update_mode = 0);

    const Eigen::Matrix2Xd& getUV() const { return UV; }
    double getEnergy() const { return total_energy; }

    void setEnergyType(EnergyType type) { energy_type = type; }
    void setPenaltyScale(double s) { penalty_scale = s; }
    void setIsoAlpha(double alpha) { iso_alpha = alpha; }
    void setLogBarrierLambda(double lambda) { log_barrier_lambda = lambda; }
    void setAMIPSVariant(AMIPSVariant v) { amips_variant = v; }

    void setLineSearchParams(double max_step, double shrink_factor) {
        ls_max_step = max_step;
        ls_shrink_factor = shrink_factor;
    }

private:
    Eigen::Matrix3Xd V;
    Eigen::Matrix3Xi F;
    mutable Eigen::Matrix2Xd UV;

    int nV = 0;
    int nF = 0;

    std::vector<Constraint> constraints;
    std::vector<bool> is_fixed;
    std::vector<Eigen::Vector2d> constraint_initial_uv;
    std::vector<Eigen::Vector2d> constraint_active_uv;
    std::vector<Eigen::Vector2d> constraint_target_uv;
    double constraint_progress = 0.0;

    struct FaceData {
        Eigen::Vector3i vid;
        Eigen::Matrix2d Dm;
        Eigen::Matrix2d Dm_inv;
        double area = 0.0;
        bool valid = true;
    };

    std::vector<FaceData> face_data;

    double total_energy = 0.0;

    EnergyType energy_type = EnergyType::AMIPS;
    AMIPSVariant amips_variant = AMIPSVariant::STANDARD;
    double penalty_scale = 10.0;
    double iso_alpha = 0.5;

    double log_barrier_lambda = 0.0;

    double ls_max_step = 1.0;
    double ls_shrink_factor = 0.85;

    std::vector<std::vector<int>> vertex_to_faces;
    std::vector<int> vertex_colors;

    double E_last_sweep = 1e20;

    bool buildLocal2DReference(int dim);
    bool initializeUVWithBoundary();
    void applyConstraints();
    void buildOneRingAdjacency();
    void computeVertexColors();

    Eigen::Matrix2d computeJacobian(int fid) const;
    double computeFaceEnergy(int fid) const;
    double computeFaceMIPS(int fid) const;
    double computeFaceAMIPS(int fid) const;
    double computeFaceIsometric(int fid) const;
    double computeAMIPSWeight(int fid) const;

    Eigen::Matrix<double, 2, 3> computeFaceGradient(int fid);
    Eigen::Matrix<double, 2, 3> computeFaceBaseGradient(int fid);
    Eigen::Matrix<double, 2, 3> numericalFaceGradient(int fid, double eps = 1e-6) const;

    bool isFaceFlipped(int fid) const;

    double computeTotalEnergy() const;
    void computeGradient(Eigen::Matrix2Xd& grad) const;

    Eigen::Vector2d computeVertexGradient(int vid) const;

    int gradientDescentSweep_GaussSeidel(double step_base);
    int gradientDescentSweep_Jacobi(double step_base);
    int gradientDescentSweep_Colored(double step_base);
    int gradientDescentSweep(double step_base);
    int gradientDescentSweep(double step_base, int mode);

    double computeVertexOneRingBoundary(int vid) const;
    bool wouldCauseFlip(int vid, const Eigen::Vector2d& new_pos) const;
    bool lineSearchForVertex(int vid, const Eigen::Vector2d& grad, double& step_size);

    void setConstraintProgress(double t);
    bool advanceConstraintProgressSoft(double target_t, double E_bound,
                                        int repair_iters = 20, int max_attempts = 10);
    bool hasAnyFlippedFace() const;
};

#endif // MIPS_H
