#include <vector>
#include <string>
#include <Eigen/Dense>
#include <limits>

class MIPS {
public:
    using Constraint = std::pair<int, Eigen::Vector3d>; 

    MIPS(const Eigen::Matrix3Xd& V,
        const Eigen::Matrix3Xi& F,
        const std::vector<Constraint>& constraints);

    ~MIPS() = default;


    // ========= 求解 =========
    bool solve(int stages,int inner_iters,int final_iters,
               double step_size = 1e-2,
               double tol = 1e-8);

    // ========= 输出 =========
    const Eigen::Matrix2Xd& getUV() const { return UV; }
    double getEnergy() const { return total_energy; }

    // ========= AMIPS 参数 =========
    enum class EnergyType { MIPS, AMIPS, ISOMETRIC };
    void setEnergyType(EnergyType type) { energy_type = type; }
    void setPenaltyScale(double s) { penalty_scale = s; }
    void setIsoAlpha(double alpha) { iso_alpha = alpha; }

private:

    Eigen::Matrix3Xd V;   // 3 x n
    Eigen::Matrix3Xi F;   // 3 x m
    Eigen::Matrix2Xd UV;  // 2 x n 当前2D参数化结果

    int nV = 0;
    int nF = 0;

    // 约束
    std::vector<Constraint> constraints;
    std::vector<bool> is_fixed;
    std::vector<Eigen::Vector2d> constraint_initial_uv;
    std::vector<Eigen::Vector2d> constraint_active_uv;
    std::vector<Eigen::Vector2d> constraint_target_uv;
    double constraint_progress = 0.0;

    // ==============================
    // 三角形预处理缓存
    // ==============================
    struct FaceData {
        Eigen::Vector3i vid;   // 三角形顶点索引
        Eigen::Matrix2d Dm;    // 参考2D局部坐标矩阵 [p1-p0, p2-p0]
        Eigen::Matrix2d Dm_inv;
        double area = 0.0;     // 参考三角形面积
        bool valid = true;
    };

    std::vector<FaceData> face_data;

    // 优化状态
    double total_energy = 0.0;

    // AMIPS 参数
    EnergyType energy_type = EnergyType::AMIPS;
    double penalty_scale = 0.0;   // 默认 0：面积加权（稳定）；>0：指数加权（压制最大失真，2D:5, 3D:2）
    double iso_alpha = 0.5;

private:
    // ========= 初始化子步骤 =========
    bool buildLocal2DReference(int dim);      // 每个三角形建立局部2D参考坐标
    bool initializeUVWithBoundary();   // 初始UV
    void applyConstraints();           // 把约束写进 UV

    // ========= 单三角形相关 =========
    Eigen::Matrix2d computeJacobian(int fid) const;
    double computeFaceEnergy(int fid) const;
    double computeFaceMIPS(int fid) const;          // 原始 MIPS
    double computeFaceAMIPS(int fid) const;          // 对称保形
    double computeFaceIsometric(int fid) const;      // isometric 能量
    double computeFaceIsoGradWeight(int fid) const;  // exp(s*E) 用于梯度加权
    Eigen::Matrix<double, 2, 3> computeFaceGradient(int fid) const;
    bool isFaceFlipped(int fid) const;

    // ========= 全局能量 =========
    double computeTotalEnergy() const;
    void computeGradient(Eigen::Matrix2Xd& grad) const;

    // ========= 优化 =========
    bool gradientDescentStep(double step_size);
    double backtrackingLineSearch(const Eigen::Matrix2Xd& grad,
                                  double init_step = 1e-2,
                                  double shrink = 0.5,
                                  int max_trials = 20);
    bool hasAnyFlippedFace() const;
    void setConstraintProgress(double t);
    bool advanceConstraintProgressSafely(double target_t,
                                         int max_backtracks = 20,
                                         double shrink = 0.5);

    // ========= 工具 =========
    Eigen::Matrix<double, 2, 3> numericalFaceGradient(int fid, double eps = 1e-6) const;
};