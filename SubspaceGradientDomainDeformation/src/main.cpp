#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include "mesh_io.h"
#include "util.h"
#include "linear-lplc.h"
#include "readcsv.h"
#include "LPLC.h"
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
        std::cerr << "can not find point" << std::endl;
        return 0;
    }
    std::vector<p_pair> pc;
    readPointCSV(argv[2], pc);
    
    vector<Vector3d> vertices;
    for (int i = 0; i < V.cols(); ++i) {
        vertices.push_back( V.col(i) );
    }
    // linearLPLC lplc(vertices, F, pc);
    // lplc.solveLinearSystem();

    
    vector<int> controller_indices1 = {36,41,191,186};
    vector<Vector3d> ab_points1 = {Vector3d(-10,0,0), Vector3d(-2,0,0)};
    vector<int> controller_indices2 = {204,209,0,5};
    vector<Vector3d> ab_points2 = {Vector3d(2,0,0), Vector3d(10,0,0)};
    boneSegment bs1(controller_indices1, ab_points1,vertices,9);
    boneSegment bs2(controller_indices2, ab_points2,vertices,9);
    vector<boneSegment>collectGamma;
    collectGamma.push_back( bs1 );
    collectGamma.push_back( bs2 );
    int rowOffset = 0;
    vector< Triplet<double> > totalTrips;
    SparseMatrix<double> Gamma_total;
    for(auto &bone : collectGamma) {
        // 将 Triplets 复制到 total
        for (int k=0; k<bone.Gamma_.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(bone.Gamma_, k); it; ++it)
                totalTrips.emplace_back(it.row() + rowOffset, it.col(), it.value());

        rowOffset += bone.r_ - 1;
    }

    Gamma_total.resize(rowOffset, V.cols());
    Gamma_total.setFromTriplets(totalTrips.begin(), totalTrips.end());
    Gamma_total.makeCompressed();
    LPLC sm(vertices, F, pc,Gamma_total);
    // sm.setBoneConstraints(controller_indices, ab_points, 9);
    sm.solve(5);
    vector<Vector3d> out_vertices = sm.deformed_;
    cout<<out_vertices.size()<<endl;
    Eigen::Matrix3Xd V_OUT(3, out_vertices.size());
    for (int i = 0; i < out_vertices.size(); ++i) {
        V_OUT.col(i) = out_vertices[i];
    }
    if (!mesh_io_helper::save_obj(argv[3], V_OUT, F)) {
        std::cerr << "save obj file error!" << std::endl;
        return 0;
    }
    cout<<"deformation done!"<<endl;




    // vector<Vector3d> out_vertices;
    // out_vertices = lplc.deformed_;
    // Eigen::Matrix3Xd V_OUT(3, V.cols());
    // for (int i = 0; i < V.cols(); ++i) {
    //     V_OUT.col(i) = out_vertices[i];
    // }
    // if (!mesh_io_helper::save_obj(argv[3], V_OUT, F)) {
    //     std::cerr << "save obj file error!" << std::endl;
    //     return 0;
    // }
    return 1;
}