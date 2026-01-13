#include "util.h"
#include <string>
#include <iostream>

using namespace std;
using namespace Eigen;

//! bb is d * 2
void bounding_box(const MatrixXd &pts, MatrixXd &bb)
{
  bb.resize(pts.rows(), 2);
  bb.col(0) = pts.rowwise().minCoeff();
  bb.col(1) = pts.rowwise().maxCoeff();
}

void merge_bounding_box(const MatrixXd &pts, MatrixXd &bb)
{
  if(bb.rows() != pts.rows() || bb.cols() < 2) {
    bounding_box(pts, bb);
  }
  else {
    MatrixXd pbb;
    bounding_box(pts, pbb);
    for(size_t i=0; i<pbb.rows(); ++i) {
      bb(i, 0) = std::min(bb(i, 0), pbb(i, 0));
      bb(i, 1) = std::max(bb(i, 1), pbb(i, 1));
    }
  }
}

double calc_res(const Eigen::MatrixXd &bb, size_t max_d)
{
  return (bb.col(1)-bb.col(0)).maxCoeff()/max_d;
}

int check_res_arg(char *res_arg)
{
  int res;
  try {
    res = stoi(res_arg);
  }
  catch(const std::invalid_argument &e) {
    cerr << res_arg << " is not a valid resolution because it's not a integer number." << endl;
    return __LINE__;
  }
  catch(const std::out_of_range &e) {
    cerr << "The resolution " << res_arg << " is out of range." << endl;
    return __LINE__;
  }
  if(res <= 0) {
    cerr << res_arg << " is not a valid resolution because it's not positive." << endl;
    return __LINE__;
  }

  return 0;
}

int check_ratio_arg(char *ratio_arg)
{
  double ratio;
  try {
    ratio = stod(ratio_arg); 
  }
  catch(const std::invalid_argument &e) {
    cerr << ratio_arg << " is not a valid ratio because it's not a number." << endl;
    return __LINE__;
  }
  catch(const std::out_of_range &e) {
    cerr << "The ratio of |E| to |M_s| " << ratio_arg << " is out of range." << endl;
    return __LINE__;
  }
}
