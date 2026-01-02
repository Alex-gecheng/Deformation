#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include "mesh_io.h"
#include "util.h"
#include "linear-lplc.h"
#include "readcsv.h"

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "number of parameter is less 3" << argc<<std::endl;
        return 0;
    }

    std::string obj_name(argv[1]);
    std::cout << obj_name ;
    Eigen::Matrix3Xd V;
    Eigen::Matrix3Xi F;
    if (!mesh_io_helper::read_obj(obj_name,V, F)) {
        std::cerr << "read obj file error!" << std::endl;
        return 0;
    }

    //read constraints;
    std::string point_name(argv[2]);
    std::ifstream file(point_name);
    if (!file.is_open()) {
        std::cerr << "can not find point��" << std::endl;
        return 0;
    }
    std::vector<p_pair> pc;
    readPointCSV(argv[2], pc);
    
    vector<Vector3d> vertices;
    for (int i = 0; i < V.cols(); ++i) {
        vertices.push_back( V.col(i) );
    }
    linearLPLC lplc(vertices, F, pc);
    lplc.solveLinearSystem();
    // lplc.solveNonLinearSystem(10);
    vector<Vector3d> out_vertices;
    out_vertices = lplc.deformed_;
    Eigen::Matrix3Xd V_OUT(3, V.cols());
    for (int i = 0; i < V.cols(); ++i) {
        V_OUT.col(i) = out_vertices[i];
    }
    if (!mesh_io_helper::save_obj(argv[3], V_OUT, F)) {
        std::cerr << "save obj file error!" << std::endl;
        return 0;
    }
    return 1;
}