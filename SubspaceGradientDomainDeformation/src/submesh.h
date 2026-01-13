#pragma once
#include <vector>
#include <map>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

class submesh {
    public:
        submesh(const vector<Vector3d>& vertices, const Matrix3Xi& faces, int controller_number);
        void meshBuilder();
        vector<Vector3d> controller_vertices_;
        Matrix3Xi controller_faces_;
    
    private:
        int controller_number_;
        vector<Vector3d> vertices_;
        Matrix3Xi faces_;
        void chooseControllers();
        void computeVertexNormals();
};  