#include "RBF.h"

#include <iostream>
#include <time.h>

using namespace std;
using namespace Eigen;

static inline void x2UAt(const VectorXd &x, size_t N, MatrixXd &U, MatrixXd &A, Vector3d &t)
{
  for(size_t d = 0; d < 3; ++d) {
    for(size_t i = 0; i < N; ++i) {
      U(d, i) = x[i*3+d];
    }
    t[d] = x[3*N+9+d];
  }
  for(size_t c = 0; c < 3; ++c) {
    for(size_t r = 0; r < 3; ++r) {
      A(r, c) = x[3*N+c*3+r];
    }
  }
}

int RBF::set_constraints(const std::vector<cons_type>& pc)
{
  P.resize(3, pc.size());
  const size_t N = P.cols();
  for (size_t i = 0; i < N; ++i)
    P.col(i) = pc[i].first;
  A.resize(3, 3);
  U.resize(3, N);
  // build and solve the eq Mx=b
  VectorXd bases;
  MatrixXd M = MatrixXd::Zero(3 * N + 12, 3 * N + 12);
  VectorXd b = VectorXd::Zero(3 * N + 12);
  for (size_t i = 0; i < N; ++i) {
    calc_bases(P.col(i), bases, type_);
    for (size_t j = 0; j < N; ++j) {
      for (size_t d = 0; d < 3; ++d) {
        M(i * 3 + d, j * 3 + d) = bases[j];
      }
    }
    for (size_t d = 0; d < 3; ++d) {
      for (size_t Ai = 0; Ai < 3; ++Ai) {
        M(i * 3 + d, N * 3 + d + Ai * 3) = P(Ai, i);
        M(N * 3 + d + Ai * 3, i * 3 + d) = P(Ai, i);
      }
    }
    for (size_t ti = 0; ti < 3; ++ti) {
      M(i * 3 + ti, N * 3 + 9 + ti) = 1;
      M(N * 3 + 9 + ti, i * 3 + ti) = 1;
    }
    for (size_t d = 0; d < 3; ++d)
      b[i * 3 + d] = pc[i].second[d];
  }

  VectorXd x;
  if (1) {
    x = M.lu().solve(b);
  }
  else {
    JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    x = svd.solve(b);
  }
  x2UAt(x, N, U, A, t);
  return 0;
}

int RBF::set_constraints(const std::vector<cons_type> &pc,
                         const std::vector<std::pair<cons_type, double> > &lc,
                         double lambda, double eta)
{
	clock_t begin = clock();
  P.resize(3, pc.size());
  const size_t N = P.cols();
  for(size_t i = 0; i < N; ++i)
    P.col(i) = pc[i].first;
  A.resize(3, 3);
  U.resize(3, N);
  // build and solve the eq Mx=b
  VectorXd bases;
  MatrixXd M = MatrixXd::Zero(3*N+12, 3*N+12);
  VectorXd b = VectorXd::Zero(3*N+12);
  for(size_t i = 0; i < N; ++i) {
    calc_bases(P.col(i), bases, type_);
    for(size_t j = 0; j < N; ++j) {
      for(size_t d = 0; d < 3; ++d) {
        M(i*3+d, j*3+d) = bases[j]-lambda*(i==j?1:0);
      }
    }
    for(size_t d = 0; d < 3; ++d) {
      for(size_t Ai = 0; Ai < 3; ++Ai) {
        M(i*3+d, N*3+d+Ai*3) = P(Ai, i);
        M(N*3+d+Ai*3, i*3+d) = P(Ai, i);
      }
    }
    for(size_t ti = 0; ti < 3; ++ti) {
      M(i*3+ti, N*3+9+ti) = 1;
      M(N*3+9+ti, i*3+ti) = 1;
    }
    for(size_t d = 0; d < 3; ++d)
      b[i*3+d] = pc[i].second[d];
  }

  VectorXd x;
  if(1) {
    x = M.lu().solve(b);
  }
  else {
    JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    x = svd.solve(b);
  }
  x2UAt(x, N, U, A, t);

  cout << "solve done: " << N << endl;
  if(lc.empty()) // no length
    return 0;

  Eigen::MatrixXd L = MatrixXd::Zero(3*lc.size(), 3*N+9+3);

  vector<vector<Eigen::VectorXd> > lc_bases(2, vector<Eigen::VectorXd>(lc.size()));
  for(size_t i = 0; i < lc.size(); ++i) { // second-first
    VectorXd phi[2];
    calc_bases(lc[i].first.first, lc_bases[0][i], type_);
    calc_bases(lc[i].first.second, lc_bases[1][i], type_);
    bases = lc_bases[1][i]-lc_bases[0][i];
    for(size_t j = 0; j < N; ++j) {
      for(size_t d = 0; d < 3; ++d)
        L(i*3+d, j*3+d) = bases[j];
    }
    Vector3d Dpi = lc[i].first.second - lc[i].first.first;
    for(size_t d = 0; d < 3; ++d) {
      for(size_t Ai = 0; Ai < 3; ++Ai) {
        L(i*3+d, N*3+Ai*3+d) = Dpi[Ai];
      }
    }
  }
  MatrixXd M2 = MatrixXd::Zero(3*N+24, 3*N+24);
  M2.block(0, 0, 3*N+12, 3*N+12) = M.block(0, 0, 3*N, 3*N+12).transpose()*M.block(0, 0, 3*N, 3*N+12);
  eta *= M2.trace()/(L.transpose()*L).trace();
  M2.block(0, 0, 3*N+12, 3*N+12) += eta*L.transpose()*L;
  M2.block(3*N+12, 0, 12, 3*N) = M.block(3*N, 0, 12, 3*N);
  M2.block(0, 3*N+12, 3*N, 12) = M.block(0, 3*N, 3*N, 12);
  const PartialPivLU<MatrixXd> lu = M2.lu();
  VectorXd b2 = VectorXd::Zero(3*N+24), x2 = VectorXd::Zero(3*N+24), b20 = M.block(0, 0, 3*N, 3*N+12).transpose()*b.head(3*N);;
  VectorXd ld(3*lc.size());
  x2.head(3*N+12) = x;
  const size_t MAX_ITR = 1000;
  size_t k = 0;
  double step = 0;
  for(; k < MAX_ITR; ++k) {
    for(size_t i = 0; i < lc.size(); ++i) {
      Vector3d ldi = L.block(i*3, 0, 3, L.cols())*x2.head(3*N+12);
      const double len = ldi.norm(), EPS = 1e-8;
      if(len > EPS)
        ldi *= lc[i].second/len;
      else {
        cout << "degenerated len" << endl;
      }
      for(size_t d = 0; d < 3; ++d)
        ld[i*3+d] = ldi[d];
    }
    b2.head(3*N+12) = b20+eta*(L.transpose()*ld);
    const VectorXd x0 = x2.head(3*N+12);
    x2 = lu.solve(b2);
    step = (x2.head(3*N+12)-x0).norm()/x0.norm();
    if(step < 1e-5)
      break;
  }
  std::cout << "cost time:" << clock() - begin <<" ms"<<std::endl;
  x2UAt(x2, N, U, A, t);
  cout << step << ":" << k << "/" << MAX_ITR << endl;
  return 0;
}

int RBF::convert(const Eigen::Vector3d &from, Eigen::Vector3d &to, double beta) const
{
  if(P.cols() == 0)
    return __LINE__;
  VectorXd bases;
  calc_bases(from,bases, type_);
  to = A*from+t+(1-beta)*U*bases;
  return 0;
}

int RBF::convert(const std::vector<Eigen::Vector3d> &from, std::vector<Eigen::Vector3d> &to, double beta) const
{
  if(P.cols() == 0)
    return __LINE__;
  to.resize(from.size());
  VectorXd bases;
  for(size_t i = 0;i < from.size();++i){
    calc_bases(from[i], bases, type_);
    to[i] = A*from[i]+t+(1-beta)*U*bases;
  }
  return 0;
}

void RBF::calc_bases(const Eigen::Vector3d &p, Eigen::VectorXd &bases) const
{
  const size_t N = P.cols();
  bases.resize(N);
  for(size_t i = 0; i < N; ++i){
    bases[i] = (P.col(i)-p).norm();
  }
}

void RBF::calc_bases(const Eigen::Vector3d &p, Eigen::VectorXd &bases, int type) const
{
  const size_t N = P.cols();
  bases.resize(N);
  switch (type) {
  case 0: // ||p_i - p||
    for(size_t i = 0; i < N; ++i){
      bases[i] = (P.col(i)-p).norm();
    }
    break;
  case 1: // ||p_i - p||^2
    for(size_t i = 0; i < N; ++i){
      bases[i] = (P.col(i)-p).squaredNorm();
    }
    break;
  case 2: // (||p_i - p||^2 + eps)^{-0.5}
    for(size_t i = 0; i < N; ++i){
      bases[i] = std::pow((P.col(i)-p).squaredNorm() + r_, 0.5);
    }
    break;
  case 3: // (||p_i - p||^2 + eps)^{-1}
    for(size_t i = 0; i < N; ++i){
      bases[i] = std::pow((P.col(i)-p).squaredNorm() + r_, -0.5);
    }
    break;
  case 4: // (||p_i - p||^2 + eps)^{-2}
    for(size_t i = 0; i < N; ++i){
      bases[i] = std::pow((P.col(i)-p).squaredNorm() + r_, -1);
    }
    break;
  case 5: // 0.5 * ||p_i - p||^2 * ln(||p_i - p||^2)
    for(size_t i = 0; i < N; ++i){
      const double d2 = (P.col(i)-p).squaredNorm();
      bases[i] = 0.5 * d2 * log(d2);
    }
    break;
  case 6: // 0.5 * ||p_i - p||^2 * ln(||p_i - p||^2)
    for(size_t i = 0; i < N; ++i){
      const double d2 = (P.col(i)-p).squaredNorm();
      bases[i] = d2 * log(sqrt(d2));
    }
    break;
  case 7: // epx( ||p_i - p||^2 / (2 * sigma))
    for(size_t i = 0; i < N; ++i){
      bases[i] = exp(-(P.col(i)-p).squaredNorm() / (2 * r_));
    }
    break;
  case 8: // epx(-||p_i - p||^2 / (2 * sigma)) * (1 + ||p_i - p||^2 / (2 * sigma));
    for(size_t i = 0; i < N; ++i){
      const double r = (P.col(i)-p).squaredNorm() / (2 * r_);
      bases[i] = exp(-r) * (1 + r);
    }
    break;
  case 9: // (||p_i - p||^2 + eps)^{-1.5}
    for(size_t i = 0; i < N; ++i){
      bases[i] = std::pow((P.col(i)-p).squaredNorm() + r_, -1.5);
    }
    break;
  case 10: // (||p_i - p||^2 + eps)^{1.5}
    for(size_t i = 0; i < N; ++i){
      bases[i] = std::pow((P.col(i)-p).squaredNorm() + r_, 1.5);
    }
    break;
  case 11:// ln(||p_i - p||^2)
    for(size_t i = 0; i < N; ++i){
      bases[i] = log((P.col(i)-p).squaredNorm());
    }
    break;
  case 12:
      for (size_t i = 0; i < N; ++i) {
          double r = (P.col(i) - p).norm();
          bases[i] = r * r * r;
      }
  default:
    break;
  }
}

// void RBF::calc_bases(const Eigen::Vector3d &p, Eigen::VectorXd &bases, double alpha) const
// {
//     const size_t N = P.cols();
//     bases.resize(N);
//     for(size_t i = 0; i < N; ++i)
//       bases[i] = std::pow((P.col(i)-p).squaredNorm() + r_, alpha);
// }
