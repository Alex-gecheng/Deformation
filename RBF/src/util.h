#ifndef HJ_EM_COMMON_H_
#define HJ_EM_COMMON_H_

#include <Eigen/Dense>

typedef std::pair<Eigen::Vector3d, Eigen::Vector3d> p_pair;

//! bb is d * 2
void bounding_box(const Eigen::MatrixXd &pts, Eigen::MatrixXd &bb);

void merge_bounding_box(const Eigen::MatrixXd &pts, Eigen::MatrixXd &bb);

//! @brief calculate resolution from bb and max_d
double calc_res(const Eigen::MatrixXd &bb, size_t max_d);

int check_res_arg(char *res_arg);

int check_ratio_arg(char *ratio_arg);

#endif
