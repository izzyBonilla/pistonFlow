#ifndef FLOW_HPP_
#define FLOW_HPP_

#include <eigen3/Eigen/Dense>

struct flowParams {
  double re;
  double ma;
  double pr;
  double gamma;
  double R;
  double uw;
  double nu;
  double L;
  double omega;
  double cp;
};

struct integParams {
  int nx;
  int ny;
  int ngx;
  int ngy;
  int nt;
  double tf;
  double dx;
  double dy;
  double dt;
};

struct flowQuant {
  Eigen::ArrayXXd rho;
  Eigen::ArrayXXd u;
  Eigen::ArrayXXd v;
  Eigen::ArrayXXd et;
};

struct Stress {
  Eigen::ArrayXXd sig11;
  Eigen::ArrayXXd sig22;
  Eigen::ArrayXXd south;
  Eigen::ArrayXXd west;
};

struct Interps {
  double east;
  double west;
  double south;
  double north;
};

Eigen::ArrayXXd rho_rhs(struct integParams integ, struct flowQuant U);
Eigen::ArrayXXd x_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S);
Eigen::ArrayXXd y_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S);
Eigen::ArrayXXd et_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S);

Eigen::ArrayXXd sig11(struct flowParams flow, struct integParams integ, struct flowQuant U);
Eigen::ArrayXXd sig22(struct flowParams flow, struct integParams integ, struct flowQuant U);
Eigen::ArrayXXd sig_south(struct flowParams flow, struct integParams integ, struct flowQuant U);
Eigen::ArrayXXd sig_west(struct flowParams flow, struct integParams integ, struct flowQuant U);

double temp(struct flowParams flow, double et, double u, double v);

#endif // FLOW_HPP_