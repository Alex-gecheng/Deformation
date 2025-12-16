#ifndef EM_REG_RBF_H_
#define EM_REG_RBF_H_

#include <vector>
#include "util.h"

class RBF
{
public:
    RBF(double r,  double alpha = 0.5, int type = 0):r_(r),alpha_(alpha), type_(type){
    }
  typedef p_pair cons_type;
  //! @notice when some points are too close, it may fail (return
  //! non-zero value) because of the numerical degenerate.
  //! 
  int set_constraints(const std::vector<cons_type>& pc);
  int set_constraints(const std::vector<cons_type> &pc,
                      const std::vector<std::pair<cons_type, double> > &lc, // length constraints
                      double lambda = 0, double eta = 0);
  //! @param beta \in [0,1], larger beta leads smoother result
  int convert(const Eigen::Vector3d &from, Eigen::Vector3d &to, double beta = 0) const;
  int convert(const std::vector<Eigen::Vector3d> &from, std::vector<Eigen::Vector3d> &to, double beta = 0) const;
private:
  void calc_bases(const Eigen::Vector3d &p, Eigen::VectorXd &bases) const;
  void calc_bases(const Eigen::Vector3d &p, Eigen::VectorXd &bases, int type) const;
  //  void calc_bases(const Eigen::Vector3d &p, Eigen::VectorXd &bases, double alpha) const;
  Eigen::MatrixXd A, U, P;
  Eigen::Vector3d t;
  int    type_;
  double alpha_, r_;
};


#endif
