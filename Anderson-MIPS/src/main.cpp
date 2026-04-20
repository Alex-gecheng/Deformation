#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include "mesh_io.h"
#include "readcsv.h"
#include "mips.h"

typedef std::pair<int, Eigen::Vector3d> p_pair;

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "number of parameter is less 4" << argc << std::endl;
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
        std::cerr << "can not find point" << std::endl;
        return 0;
    }
    std::vector<p_pair> pc;
    readPointCSV(argv[2], pc);
    MIPS solver(V, F, pc);


    if (!solver.solve(200,15,100, 1e-2,1e-6)) {
        std::cerr << "MIPS solve failed.\n";
        return 0;
    }

    const Eigen::Matrix2Xd& UV = solver.getUV();

    // TODO: 保存 UV 到 obj / txt
    std::cout << "Done. UV cols = " << UV.cols() << std::endl;

    

    if (!mesh_io_helper::save_obj(argv[3], UV, F)) {
        std::cerr << "save obj file error!" << std::endl;
        return 0;
    }
    return 1;
}