#ifndef MESH_IO_H_H_
#define MESH_IO_H_H_

#include <string>
#include <Eigen/Dense>

namespace  mesh_io_helper {
	template<class Vector, class IS>
	void read_vector_binary(IS &in, Vector& vec) {
		typename Vector::Index size = 0;
		in.read((char*)(&size), sizeof(typename Vector::Index));
		vec.resize(size);
		in.read((char *)vec.data(), size * sizeof(typename Vector::Scalar));
	}

	template<class Matrix, class IS>
	void read_matrix_binary(IS &in, Matrix& matrix) {
		typename Matrix::Index rows = 0, cols = 0;
		in.read((char*)(&rows), sizeof(typename Matrix::Index));
		in.read((char*)(&cols), sizeof(typename Matrix::Index));
		matrix.resize(rows, cols);
		in.read((char*)matrix.data(), rows*cols * sizeof(typename Matrix::Scalar));
	}

	
	// @brief read and write obj mesh
// @param filename: the file name of the mesh
// @param V: the vertex list of the mesh
// @param F: the face list of the mesh
// @param tV: the texture vertex list of the mesh
// @param tF: the texture face list of the mesh
// @return 0: read or write failed; 1: read or write successufully

int read_obj(const std::string &filename, Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F);
int save_obj(const std::string &filename, const Eigen::MatrixXd &nods,const Eigen::Matrix3Xi &tris);

int read_obj(const std::string &in_file,
             Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F,
             Eigen::Matrix2Xd &tV, Eigen::Matrix3Xi &tF);
int save_obj(const std::string &out_file,
             const Eigen::Matrix3Xd &V,const Eigen::Matrix3Xi &F,
             const Eigen::Matrix2Xd &tv,const Eigen::Matrix3Xi &tf);

int read_off(const std::string &filename,
	         Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F);

int read_tri_from_vtk(const std::string &path,Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F);				 
int read_tet_from_vtk(const std::string &path,Eigen::Matrix3Xd &V, Eigen::Matrix4Xi &T);

}


#endif
