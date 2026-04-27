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
    bool solve(int stages,int inner_iters,
               double step_size = 1e-2,
               double tol = 1e-8);

    // ========= 输出 =========
    const Eigen::Matrix2Xd& getUV() const { return UV; }
    double getEnergy() const { return total_energy; }

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

private:
    // ========= 初始化子步骤 =========
    bool buildLocal2DReference(int dim);      // 每个三角形建立局部2D参考坐标
    bool initializeUVWithBoundary();   // 初始UV
    void applyConstraints();           // 把约束写进 UV

    // ========= 单三角形相关 =========
    Eigen::Matrix2d computeJacobian(int fid) const;
    double computeFaceEnergy(int fid) const;
    Eigen::Matrix<double, 2, 3> computeFaceGradient(int fid) const;
    bool isFaceFlipped(int fid) const;

    // ========= 能量 =========
    double computeTotalEnergy() const;
    Eigen::Vector2d computeLocalGradient(int v) const;
    double computeLocalEnergy(int v) const;
    
    // ========= 优化 =========
    bool hasAnyFlippedFace() const;
    void setConstraintProgress(double t);
    bool advanceConstraintProgressSafely(double target_t,
                                         int max_backtracks = 20,
                                         double shrink = 0.5);
    void update_vertex_gd(int v, double init_step = 1e-2);

    // ========= 工具 =========
    Eigen::Matrix<double, 2, 3> numericalFaceGradient(int fid, double eps = 1e-6) const;

private:
    // ========= 图着色 =========
    // 图结构
    std::vector<std::vector<int>> vertex_adj;
    std::vector<std::vector<int>> vertex_faces;

    // 图着色
    std::vector<int> vertex_color;
    std::vector<std::vector<int>> color_groups;
    int num_colors = 0;
    void buildVertexAdjacency();
    void buildVertexFaceAdjacency();
    void buildGraphColoring();
};