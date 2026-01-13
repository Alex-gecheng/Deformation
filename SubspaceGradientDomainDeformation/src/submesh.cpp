#include "submesh.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;


submesh::submesh(const vector<Vector3d>& vertices, const Matrix3Xi& faces, int controller_number = 300)
    : vertices_(vertices), faces_(faces), controller_number_(controller_number) {
    meshBuilder();
}

void submesh::meshBuilder() {
    chooseControllers();
    // computeVertexNormals();

}

// void submesh::chooseControllers() {
//     controller_vertices_.clear();
//     //采样密度
//     int space = vertices_.size() / controller_number_;
//     for (int i = 0; i < vertices_.size(); i += space) {
//         controller_vertices_.push_back(vertices_[i]);
//         if (controller_vertices_.size() >= controller_number_) break;
//     }
//     // 构建controller_faces_
//     std::vector<Point_3> points;
//     for (const auto& v : controller_vertices_) {
//         points.push_back(Point_3(v.x(), v.y(), v.z()));
//     }

//     Polyhedron poly;
//     CGAL::convex_hull_3(points.begin(), points.end(), poly);

//     // 三角面数
//     int n_faces = std::distance(poly.facets_begin(), poly.facets_end());
//     controller_faces_.resize(3, n_faces);

//     // 填充 faces
//     int f_idx = 0;
//     for (auto f = poly.facets_begin(); f != poly.facets_end(); ++f, ++f_idx) {
//         auto h = f->facet_begin();
//         int v0 = std::distance(poly.vertices_begin(), h->vertex()); ++h;
//         int v1 = std::distance(poly.vertices_begin(), h->vertex()); ++h;
//         int v2 = std::distance(poly.vertices_begin(), h->vertex());
//         controller_faces_.col(f_idx) << v0, v1, v2;
//     }
// }

void submesh::chooseControllers()
{
    controller_vertices_.clear();
    controller_faces_.resize(3, 0);

    // =========================
    // 1. 均匀下采样顶点
    // =========================
    int n = vertices_.size();
    if (n == 0 || controller_number_ <= 0) return;

    int step = std::max(1, n / controller_number_);
    for (int i = 0; i < n && controller_vertices_.size() < controller_number_; i += step) {
        controller_vertices_.push_back(vertices_[i]);
    }

    // =========================
    // 2. 构造 CGAL 点集
    // =========================
    std::vector<Point_3> points;
    points.reserve(controller_vertices_.size());
    for (const auto& v : controller_vertices_) {
        points.emplace_back(v.x(), v.y(), v.z());
    }

    // =========================
    // 3. 计算凸包
    // =========================
    Polyhedron poly;
    CGAL::convex_hull_3(points.begin(), points.end(), poly);

    if (poly.empty()) return;

    // =========================
    // 4. 建立 poly vertex -> controller index 映射
    // =========================
    std::unordered_map<
        Polyhedron::Vertex_const_handle,
        int
    > poly_to_ctrl;

    // ⚠️ 凸包会丢点，这里只保留 hull 上的点
    int idx = 0;
    controller_vertices_.clear();  // 只保留 hull vertices

    for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v, ++idx) {
        const Point_3& p = v->point();
        controller_vertices_.emplace_back(
            p.x(), p.y(), p.z()
        );
        poly_to_ctrl[v] = idx;
    }

    // =========================
    // 5. 统计三角面数量（facet 可能是 n-gon）
    // =========================
    int tri_count = 0;
    for (auto f = poly.facets_begin(); f != poly.facets_end(); ++f) {
        int deg = 0;
        auto h = f->facet_begin();
        do {
            ++deg;
            ++h;
        } while (h != f->facet_begin());

        if (deg >= 3) tri_count += deg - 2;
    }

    controller_faces_.resize(3, tri_count);

    // =========================
    // 6. 面 triangulation（fan）
    // =========================
    int f_idx = 0;
    for (auto f = poly.facets_begin(); f != poly.facets_end(); ++f) {

        // 收集一个 facet 的所有顶点
        std::vector<int> face_indices;
        auto h = f->facet_begin();
        do {
            face_indices.push_back(poly_to_ctrl[h->vertex()]);
            ++h;
        } while (h != f->facet_begin());

        // fan triangulation
        for (int i = 1; i + 1 < face_indices.size(); ++i) {
            controller_faces_.col(f_idx++) <<
                face_indices[0],
                face_indices[i],
                face_indices[i + 1];
        }
    }
}
