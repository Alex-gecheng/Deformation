#include "mesh_io.h"
#include <fstream>
#include <iostream>
#include <vector>

namespace mesh_io_helper {
int read_obj(const std::string &filename,
            Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F)
{
    std::ifstream is(filename.c_str());
	if (!is)
		return 0;  
    std::vector<double>  vs;
    std::vector<int>  fs;
    std::string line, pair[3];
    double  node[3];
    int  tri;
    while (!is.eof()) {
        std::getline(is, line);
        if (line.empty() || 13 == line[0])
            continue;
        std::istringstream instream(line);
        std::string word;
        instream >> word;
        if ("v" == word || "V" == word) {
            instream >> node[0] >> node[1] >> node[2];
            for (size_t j = 0; j < 3; ++j){
                vs.push_back(node[j]);
            }
        }
        else if ('f' == word[0] || 'F' == word[0]) {
            instream >> pair[0] >> pair[1] >> pair[2];
            for (size_t j = 0; j < 3; ++j) {
                tri = strtoul(pair[j].c_str(),NULL,10) - 1;
                fs.push_back(tri);
            }
        }
    }
    is.close();
	V = Eigen::Map<Eigen::Matrix3Xd>(vs.data(), 3, vs.size() / 3);
	F = Eigen::Map<Eigen::Matrix3Xi>(fs.data(), 3, fs.size() / 3);
    return 1;
}

int save_obj(const std::string &filename,
            const Eigen::MatrixXd &nods,
            const Eigen::Matrix3Xi &tris)
{
    std::ofstream os(filename.c_str());
	if (!os)
		return 0;
	if (nods.rows() == 2) {
		for (size_t i = 0; i < nods.cols(); ++i) {
			os << "v " << nods(0, i) <<" "<<nods(1, i)<<" 0\n";
		}
	}
	else {
		for (size_t i = 0; i < nods.cols(); ++i) {
			os << "v " << nods.col(i).transpose() << "\n";
		}
	}
    for (size_t i = 0; i < tris.cols(); ++i){
        const Eigen::Vector3i f = tris.col(i) + Eigen::Vector3i::Ones();
        os << "f " << f.transpose() << "\n";
    }
    os.close();
    return 1;
}

int read_obj(const std::string &in_file,
             Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F,
             Eigen::Matrix2Xd &tV, Eigen::Matrix3Xi &tF)
{
    std::ifstream is(in_file);
    if(!is)
        return 0;
    std::vector<double>  vs, vt;
    std::vector<int>  fs, ft;
    std::string line, pair[3];
    double  node[3];
    while (!is.eof()) {
        std::getline(is, line);
        if (line.empty() || 13 == line[0])
            continue;
        std::istringstream instream(line);
        std::string word;
        instream >> word;
        if ("v" == word || "V" == word) {
            instream >> node[0] >> node[1] >> node[2];
            for (size_t j = 0; j < 3; ++j){
                vs.push_back(node[j]);
            }
        }else if ("vt" == word || "VT" == word) {
            instream >> node[0] >> node[1];
            for (size_t j = 0; j < 2; ++j){
                vt.push_back(node[j]);
            }
        }
        else if ('f' == word[0] || 'F' == word[0]) {
            instream >> pair[0] >> pair[1] >> pair[2];
            //get vertex id in the this triangle face
            size_t tri;
            for (size_t j = 0; j < 3; ++j) {
                std::string vertex_str = pair[j].substr(0, pair[j].find('/'));
                sscanf(vertex_str.c_str(), "%lu", &tri);
                fs.push_back(--tri);

                vertex_str = pair[j].substr(pair[j].find('/')+1);
                sscanf(vertex_str.c_str(), "%lu", &tri);
                ft.push_back(--tri);
            }
        }
    }
	V = Eigen::Map<Eigen::Matrix3Xd>(vs.data(), 3, vs.size() / 3);
	F = Eigen::Map<Eigen::Matrix3Xi>(fs.data(), 3, fs.size() / 3);
	tV = Eigen::Map<Eigen::Matrix2Xd>(vt.data(), 2, vt.size() / 2);
	tF = Eigen::Map<Eigen::Matrix3Xi>(ft.data(), 3, ft.size() / 3);
    return 1;
}

int save_obj(const std::string &out_file,
             const Eigen::Matrix3Xd &V,const Eigen::Matrix3Xi &F,
             const Eigen::Matrix2Xd &tv,const Eigen::Matrix3Xi &tf)
{
    std::ofstream os(out_file);
    if(!os)
        return 0;
    os <<"mtllib ./make_human.mtl"<<std::endl;
    for (int i = 0; i < V.cols(); ++i){
        os << "v " << V(0,i) << " " << V(1,i) << " " << V(2,i) << "\n";
    }
    for (size_t i = 0; i < tv.cols(); ++i){
        os << "vt " << tv(0, i) << " " << tv(1, i)<< "\n";
    }
    for (size_t i = 0; i < tf.cols(); ++i){
        os << "f " << F(0, i) + 1 << "/" << tf(0, i) + 1 << " "
           << F(1, i)  + 1 << "/" << tf(1, i) + 1 << " "
           << F(2, i)  + 1 << "/" << tf(2, i) + 1 << "\n";
    }
    os.close();
    return 1;
}

int read_off(const std::string &filename,
	         Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F)
{
	std::ifstream infile(filename);
	if (!infile)
		return 0;
	std::string temp;
	infile >> temp;
	size_t numVertices, numFaces;
	infile >> numVertices >> numFaces >> temp;

	if (infile.eof())
		return 0;
	V.resize(3, numVertices);
	for (size_t i = 0; i < numVertices; i++) {
		if (infile.eof())
			return 0;
		infile >> V(0, i) >> V(1, i) >> V(2, i);
	}

	int three;
	F.resize(3, numFaces);
	for (unsigned int i = 0; i < numFaces; i++) {
		if (infile.eof()) return 0;
		infile >> three >> F(0, i) >> F(1, i) >> F(2, i);
	}
	return 1;
}

int read_tri_from_vtk(const std::string &path,Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F)
{
    std::ifstream ifs(path);
    if(ifs.fail()) {
        std::cerr << "[info] " << "can not open file" << path << std::endl;
        return __LINE__;
    }

    std::string str;
    size_t point_num = 0,cell_num = 0;
    while(!ifs.eof()){
        ifs >> str;
        if(str == "POINTS"){
            ifs >> point_num >> str;
            V.resize(3, point_num);
            for(size_t i = 0;i < point_num; ++i){
                for(size_t j = 0;j < 3; ++j)
                    ifs >> V(j, i);
            }
            continue;
        }
        if(str == "CELLS"){
            ifs >> cell_num >> str;
            size_t point_number_of_cell = 0;
			F.resize(3, cell_num);
            for(size_t ci = 0; ci < cell_num; ++ci){
                ifs >> point_number_of_cell;
                if(point_number_of_cell != 3){
                    for(size_t i = 0; i < point_number_of_cell; ++i)
                        ifs >> str;
                }else{
                    for(size_t i = 0; i < point_number_of_cell; ++i){
						ifs >> F(i, ci);
                    }
                }
            }
        }
    }
    return 0;
}

int read_tet_from_vtk(const std::string &path,Eigen::Matrix3Xd &V, Eigen::Matrix4Xi &T)
{
    std::ifstream ifs(path);
    if(ifs.fail()) {
        std::cerr << "[info] " << "can not open file" << path << std::endl;
        return __LINE__;
    }

    std::string str;
    size_t point_num = 0,cell_num = 0;
    while(!ifs.eof()){
        ifs >> str;
        if(str == "POINTS"){
            ifs >> point_num >> str;
            V.resize(3, point_num);
            for(size_t i = 0;i < point_num; ++i){
                for(size_t j = 0;j < 3; ++j)
                    ifs >> V(j, i);
            }
            continue;
        }
        if(str == "CELLS"){
            ifs >> cell_num >> str;
            size_t point_number_of_cell = 0;
			T.resize(4, cell_num);
            for(size_t ci = 0; ci < cell_num; ++ci){
                ifs >> point_number_of_cell;
                if(point_number_of_cell != 4){
                    for(size_t i = 0; i < point_number_of_cell; ++i)
                        ifs >> str;
                }else{
                    for(size_t i = 0; i < point_number_of_cell; ++i){
                        ifs >> T(i, ci);
                    }
                }
            }
        }
    }
    return 0;
}
}
