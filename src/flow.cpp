#include <fstream>
#include <iostream>
#include <cmath>

#include "flow.hpp"

#define TI 300 // 300 Kelvin initial temperature
#define PI 100000 // 100000 Pa
#define NGHOST 2
#define MA0 0.025 // default mach

using Eigen::ArrayXXd;

int main(int argc, char* argv[]) {

  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n");

  // define structs to carry important quantitites
  struct flowParams n; // parameters defining flow like Re and Ma
  struct integParams integ; // integrator parameters
  struct flowQuant U; // flow quantities: density, u velocity, v velocity, total energy
  struct Stress S; // stress grids

  if(argc == 2) {
    n.ma = atof(argv[1]);
  } else if (argc == 1) {
    n.ma = MA0;
  } else {
    std::cout << "Please either input mach number or no arguments for Ma = 0.025";
    exit(1);
  }


  // input flow quantities
  n.re = 100;
  n.pr = 0.7;
  n.gamma = 1.4;
  n.R = 287;
  double a = sqrt(n.gamma*n.R*TI);
  n.nu = MA0*a/n.re;
  n.uw = n.ma*a;
  n.L = n.re*n.nu/n.uw;
  n.nu = n.uw*n.L/n.re;
  n.cp = n.gamma/(n.gamma-1)*n.R;
  n.omega = pow(1/n.L,2)*2*n.nu; // * default value of this is 0.173594

  // input integrator quantities
  integ.tf = std::round(10/n.omega*100.0)/100.0;
  integ.dt = 0.0000001;
  integ.nt = integ.tf/integ.dt;

  std::cout << "tf = " << integ.tf << std::endl;
  std::cout << "L = " << n.L << std::endl;
  std::cout << "Uw = " << n.uw << std::endl;

  integ.nx = 40;
  integ.ngx = integ.nx + NGHOST;
  integ.ny = 40;
  integ.ngy = integ.ny + NGHOST;

  integ.dx = n.L/integ.nx;
  integ.dy = n.L/integ.ny;

  // get initial density and enery of flow field
  double rho_i = PI/(n.R*TI);
  double et_i = n.R*TI/(n.gamma-1);

  // initialize flow quantity arrays
  U.rho = ArrayXXd::Constant(integ.ngx,integ.ngy,rho_i);
  U.et = ArrayXXd::Constant(integ.ngx,integ.ngy,et_i);

  U.u = ArrayXXd::Zero(integ.ngx,integ.ngy); // * flow is quiescent to start
  U.v = ArrayXXd::Zero(integ.ngx,integ.ngy);

  // ! Array indices indicate direction such that increasing k represents increasing
  // ! x while increasing l represents increasing y


  // ! Stress scheme is to store the -1/2th stress at k,l
  // ! i.e. sig_11|k-1/2,l is stored at S.sig11(k,l)
  // initialize stress arrays
  S.sig11 = ArrayXXd::Zero(integ.ngx,integ.ngy);
  S.sig22 = ArrayXXd::Zero(integ.ngx,integ.ngy);

  // non-principal stresses on southern and western cell faces
  S.south = ArrayXXd::Zero(integ.ngx,integ.ngy);
  S.west = ArrayXXd::Zero(integ.ngx,integ.ngy);


  // * placeholder
  ArrayXXd rho_old = ArrayXXd::Zero(integ.ngx,integ.ngy);

  double t;
  double T_minus;
  double wall_velo;

  // time integration
  for(int s = 0; s < integ.nt; ++s) {

    t = s*integ.dt;
    wall_velo = n.uw*sin(n.omega*t);

    // handle boundary conditions first
    // leave the corner values out of this
    // * top and bottom walls
    for(int k = 1; k < integ.ngx-1; ++k) {
      // density
      U.rho(k,integ.ngy-1) = U.rho(k,integ.ngy-2);
      U.rho(k,0) = U.rho(k,1);

      // top velocity boundary conditions
      U.u(k,integ.ngy-1) = 2*wall_velo-U.u(k,integ.ngy-2);
      U.v(k,integ.ngy-1) = -U.v(k,integ.ngy-2);

      // bottom velocity boundary conditions
      U.u(k,0) = -U.u(k,1);
      U.v(k,0) = -U.v(k,1);

      // energy boundary conditions
      U.et(k,integ.ngy-1) = (2*TI-temp(n,U.et(k,integ.ngy-2),U.u(k,integ.ngy-2),U.v(k,integ.ngy-2)))*n.R/(n.gamma-1) +
        0.5*(pow(U.u(k,integ.ngy-1),2)+pow(U.v(k,integ.ngy-1),2));
      U.et(k,0) = (2*TI-temp(n,U.et(k,1),U.u(k,1),U.v(k,1)))*n.R/(n.gamma-1)+0.5*(pow(U.u(k,0),2)+pow(U.v(k,0),2)); // bottom
    }
    
    // * left and right walls
    for(int l = 1; l < integ.ngy-1; ++l) {

      // density
      U.rho(integ.ngx-1,l) = U.rho(integ.ngx-2,l);
      U.rho(0,l) = U.rho(1,l);

      // right velocity boundary conditions
      U.u(integ.ngx-1,l) = -U.u(integ.ngx-2,l);
      U.v(integ.ngx-1,l) = -U.v(integ.ngx-2,l);

      // left velocity boundary conditions
      U.u(0,l) = -U.u(1,l);
      U.v(0,l) = -U.v(1,l);

      // energy
      U.et(integ.ngx-1,l) = n.R/(n.gamma-1)*(2*TI-temp(n,U.et(integ.ngx-2,l),U.u(integ.ngx-2,l),U.v(integ.ngx-2,l))) +
        0.5*(pow(U.u(integ.ngx-1,l),2)+pow(U.v(integ.ngx-1,l),2));
      U.et(0,l) = n.R/(n.gamma-1)*(2*TI-temp(n,U.et(1,l),U.u(1,l),U.v(1,l))) +
        0.5*(pow(U.u(0,l),2)+pow(U.v(0,l),2));
    }

    S.sig11 = sig11(n,integ,U);
    S.sig22 = sig22(n,integ,U);
    S.south = sig_south(n,integ,U);
    S.west = sig_west(n,integ,U);

    rho_old = U.rho;

    U.rho = U.rho + integ.dt*rho_rhs(integ,U);
    U.u = (U.u*rho_old + integ.dt*x_rhs(n,integ,U,S))/U.rho;
    U.v = (U.v*rho_old + integ.dt*y_rhs(n,integ,U,S))/U.rho;
    U.et = (U.et*rho_old + integ.dt*et_rhs(n,integ,U,S))/U.rho;

  }

  // Temperature post-process

  U.et = (n.gamma-1)/n.R*(U.et-0.5*(U.u*U.u+U.v*U.v));

    // write density to file
  std::ofstream rho_file("rho.csv");
  if(rho_file.is_open()) {
      rho_file << U.rho.format(CSVFormat);
  } else {
      std::cout << "Can't write to file rho.csv \n";
      exit(1);
  }

  // write u to file
  std::ofstream u_file("u.csv");
  if(u_file.is_open()) {
      u_file << U.u.format(CSVFormat);
  } else {
      std::cout << "Can't write to file u.csv \n";
      exit(1);
  }

  // write v to file
  std::ofstream v_file("v.csv");
  if(v_file.is_open()) {
      v_file << U.v.format(CSVFormat);
  } else {
      std::cout << "Can't write to file v.csv \n";
      exit(1);
  }

  // write et to file
  std::ofstream et_file("et.csv");
  if(et_file.is_open()) {
      et_file << U.et.format(CSVFormat);
  } else {
      std::cout << "Can't write to file et.csv \n";
      exit(1);
  }

  return 0;
}


Eigen::ArrayXXd rho_rhs(struct integParams integ, struct flowQuant U) {
  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  int k,l;

  #pragma omp parallel for private(k) collapse(2)
  for(l = 1; l < integ.ngx-1; ++l) {
    for(k = 1; k < integ.ngx-1; ++k) {
      f(k,l) =
        -(U.rho(k+1,l)*U.u(k+1,l)-U.rho(k-1,l)*U.u(k-1,l))/(2*integ.dx)
        -(U.rho(k,l+1)*U.v(k,l+1)-U.rho(k,l-1)*U.v(k,l-1))/(2*integ.dy);
    }
  }

  return f;
}

Eigen::ArrayXXd x_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S) {
  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  int k,l;

  #pragma omp parallel for private(k) collapse(2)
  for(l = 1; l < integ.ngy-1; ++l) {
    for(k = 1; k < integ.ngx-1; ++k) {
      f(k,l) =
        - (U.rho(k+1,l)*pow(U.u(k+1,l),2)-U.rho(k-1,l)*pow(U.u(k-1,l),2))/(2*integ.dx)
        - (U.rho(k,l+1)*U.u(k,l+1)*U.v(k,l+1)-U.rho(k,l-1)*U.u(k,l-1)*U.v(k,l-1))/(2*integ.dy)
        + (S.sig11(k+1,l)-S.sig11(k,l))/integ.dx
        + (S.south(k,l+1)-S.south(k,l))/integ.dy;
    }
  }

  return f;
}

Eigen::ArrayXXd y_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S) {
  
  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  int k,l;

  #pragma omp parallel for private(k) collapse(2)
  for(l = 1; l < integ.ngy-1; ++l) {
    for(k = 1; k < integ.ngx-1; ++k) {
      f(k,l) =
        - (U.rho(k,l+1)*pow(U.v(k,l+1),2)-U.rho(k,l-1)*pow(U.v(k,l-1),2))/(2*integ.dy)
        - (U.rho(k+1,l)*U.u(k+1,l)*U.v(k+1,l)-U.rho(k-1,l)*U.u(k-1,l)*U.v(k-1,l))/(2*integ.dx)
        + (S.west(k+1,l)-S.west(k,l))/integ.dx
        + (S.sig22(k,l+1)-S.sig22(k,l))/integ.dy;
    }
  }

  return f;
}

Eigen::ArrayXXd et_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S) {
  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  // variables for holding interpolated lambda
  double l_up;
  double l_down;
  double l_left;
  double l_right;

  // variables for holding interpolated u
  double u_up;
  double u_down;
  double u_left;
  double u_right;

  // variables for holding interpolated v
  double v_up;
  double v_down;
  double v_left;
  double v_right;

  // variables for holding temperature values
  double Tc;
  double Tkp;
  double Tkm;
  double Tlp;
  double Tlm;

  int k,l;

  #pragma omp parallel for private(k,l_up,l_down,l_left, l_right, u_up, u_down, u_left) \
  private(u_right, v_up, v_down, v_left, v_right, Tc, Tkp, Tkm, Tlp, Tlm) collapse(2)
  for(l = 1; l < integ.ngy-1; ++l) {
    for(k = 1; k < integ.ngx-1; ++k) {

      // * calculate lambdas as l = rho*cp*nu/Pr
      l_up = (U.rho(k,l+1)+U.rho(k,l))/2*flow.cp*flow.nu/flow.pr;
      l_down = (U.rho(k,l-1)+U.rho(k,l))/2*flow.cp*flow.nu/flow.pr;
      l_left = (U.rho(k-1,l)+U.rho(k,l))/2*flow.cp*flow.nu/flow.pr;
      l_right = (U.rho(k+1,l)+U.rho(k,l))/2*flow.cp*flow.nu/flow.pr;

      // * interpolate velocities
      v_up = (U.v(k,l+1)+U.v(k,l))/2;
      v_down = (U.v(k,l-1)+U.v(k,l))/2;
      v_left = (U.v(k-1,l)+U.v(k,l))/2;
      v_right = (U.v(k+1,l)+U.v(k,l))/2;
  
      u_up = (U.u(k,l+1)+U.u(k,l))/2;
      u_down = (U.u(k,l-1)+U.u(k,l))/2;
      u_left = (U.u(k-1,l)+U.u(k,l))/2;
      u_right = (U.u(k+1,l)+U.u(k,l))/2;

      // * Calculate Temperatures
      Tc = temp(flow,U.et(k,l),U.u(k,l),U.v(k,l));
      Tkp = temp(flow,U.et(k+1,l),U.u(k+1,l),U.v(k+1,l));
      Tkm = temp(flow,U.et(k-1,l),U.u(k-1,l),U.v(k-1,l));
      Tlp = temp(flow,U.et(k,l+1),U.u(k,l+1),U.v(k,l+1));
      Tlm = temp(flow,U.et(k,l-1),U.u(k,l-1),U.v(k,l-1));

      f(k,l) =
        // * diffusion term
        - (U.rho(k+1,l)*U.u(k+1,l)*U.et(k+1,l) - (U.rho(k-1,l)*U.u(k-1,l)*U.et(k-1,l)))/(2*integ.dx)
        - (U.rho(k,l+1)*U.v(k,l+1)*U.et(k,l+1) - (U.rho(k,l-1)*U.v(k,l-1)*U.et(k,l-1)))/(2*integ.dx)

        // * Thermal conductivity
        + (l_right*(Tkp-Tc)+l_left*(Tkm-Tc))/(integ.dx*integ.dx)
        + (l_up*(Tlp-Tc)+l_down*(Tlm-Tc))/(integ.dy*integ.dy)

        // * Stresses
        + (S.sig11(k+1,l)*u_right+S.west(k+1,l)*v_right-S.sig11(k,l)*u_left-S.west(k,l)*v_left)/integ.dx
        + (S.sig22(k,l+1)*v_up+S.south(k,l+1)*u_up - S.sig22(k,l)*v_down - S.south(k,l)*u_down)/integ.dy;
    }
  }

  return f;
}


Eigen::ArrayXXd sig11(struct flowParams flow, struct integParams integ, struct flowQuant U) {

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  double interp_rho;
  double pressure;
  double rt;
  double mag_v;
  double mu;

  double interp_v_plus;
  double interp_v_minus;

  int k,l;

  #pragma omp parallel for private(k,interp_rho,pressure,rt,mag_v,mu,interp_v_plus,interp_v_minus) collapse(2)
  for(l = 1; l < integ.ngy-1; ++l) {
    for(k = 1; k < integ.ngx; ++k) {

      interp_rho = (U.rho(k-1,l)+U.rho(k,l))/2; // interpolate rho between k-1 and k to get k-1/2
      mag_v = pow((U.u(k-1,l)+U.u(k,l))/2,2)+pow((U.v(k-1,l)+U.v(k,l))/2,2); // interp u and v
      rt = ((U.et(k-1,l)+U.et(k,l))/2-mag_v/2)*(flow.gamma-1); // RT = (et - 0.5*|V|^2)*(gamma-1) 
      pressure = interp_rho*rt; // eos

      mu = interp_rho*flow.nu;
      interp_v_plus = (U.v(k,l)+U.v(k-1,l)+U.v(k-1,l+1)+U.v(k,l+1))/4; // top left corner v velo
      interp_v_minus = (U.v(k,l)+U.v(k,l-1)+U.v(k-1,l-1)+U.v(k-1,l))/4; // bottom left corner v velo

      sigma(k,l) = -pressure + mu * (
        4/3*(U.u(k,l)-U.u(k-1,l))/integ.dx - 2/3*(interp_v_plus-interp_v_minus)/integ.dy
      );
    }
  }

  // ! handle top corners separately
  // *  top left 
  interp_rho = (U.rho(0,integ.ngy-2)+U.rho(1,integ.ngy-2))/2;
  mu = interp_rho*flow.nu;
  sigma(1,integ.ngy-2) = -interp_rho*flow.R*TI+mu*(
    4/3*(U.u(1,integ.ngy-2)-U.u(0,integ.ngy-2))/integ.dx
  );

  // * top right
  interp_rho = (U.rho(integ.ngx-2,integ.ngy-2)+U.rho(integ.ngx-1,integ.ngy-2))/2;
  mu = interp_rho*flow.nu;
  sigma(integ.ngx-1,integ.ngy-2) = -interp_rho*flow.R*TI+mu*(
    4/3*(U.u(integ.ngx-1,integ.ngy-2)-U.u(integ.ngx-2,integ.ngy-2))/integ.dx
  );

  return sigma;
}

Eigen::ArrayXXd sig22(struct flowParams flow, struct integParams integ, struct flowQuant U) {

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  double interp_rho;
  double pressure;
  double rt;
  double mag_v;
  double mu;

  double interp_u_plus;
  double interp_u_minus;

  int k,l;

  #pragma omp parallel for private(k,interp_rho,pressure,rt,mag_v,mu,interp_u_plus,interp_u_minus) collapse(2)
  for(l = 1; l < integ.ngy; ++l) {
    for(k = 1; k < integ.ngx-1; ++k) {
      interp_rho = (U.rho(k,l)+U.rho(k,l-1))/2;
      mu = interp_rho*flow.nu;
      mag_v = pow((U.u(k,l-1)+U.u(k,l))/2,2)+pow((U.v(k,l-1)+U.v(k,l))/2,2); // interp u and v
      rt = ((U.et(k,l-1)+U.et(k,l))/2-mag_v/2)*(flow.gamma-1);

      pressure = interp_rho*rt;

      interp_u_plus = (U.u(k,l)+U.u(k+1,l)+U.u(k+1,l-1)+U.u(k,l-1))/4;
      interp_u_minus = (U.u(k,l)+U.u(k,l-1)+U.u(k-1,l-1)+U.u(k-1,l))/4;

      sigma(k,l) = -pressure + mu * (
        4/3*(U.v(k,l)-U.v(k,l-1))/integ.dy - 2/3*(interp_u_plus-interp_u_minus)/integ.dx
      );

    }
  }

  // * top left corner
  interp_rho = (U.rho(1,integ.ngy-1)+U.rho(1,integ.ngy-2))/2;
  mu = interp_rho*flow.nu;
  sigma(1,integ.ngy-1) = -interp_rho*flow.R*TI+mu*(
    4/3*(U.v(1,integ.ngy-1)-U.v(1,integ.ngy-2))/integ.dy
  );

  // * top right corner

  interp_rho = (U.rho(integ.ngx-2,integ.ngy-1)+U.rho(integ.ngx-2,integ.ngy-2))/2;
  sigma(integ.ngx-2,integ.ngy-1) = -interp_rho*flow.R*TI+mu*(
    4/3*(U.v(integ.ngx-2,integ.ngy-1)-U.v(integ.ngx-2,integ.ngy-2))/integ.dy
  );

  return sigma;
}

Eigen::ArrayXXd sig_south(struct flowParams flow, struct integParams integ, struct flowQuant U) {

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  double interp_rho;
  double mu;

  double interp_v_plus;
  double interp_v_minus;

  int k,l;

  #pragma omp for private(k,interp_rho,mu,interp_v_plus,interp_v_minus) collapse(2)
  for(l = 1; l < integ.ngy; ++l) {
    for(k = 1; k < integ.ngx-1; ++k) {
      interp_rho = (U.rho(k,l-1)+U.rho(k,l))/2;
      mu = interp_rho*flow.nu;

      interp_v_plus = (U.v(k,l)+U.v(k,l-1)+U.v(k-1,l-1)+U.v(k-1,l))/4;
      interp_v_minus = (U.v(k,l)+U.v(k+1,l)+U.v(k+1,l-1)+U.v(k,l-1))/4;

      sigma(k,l) = mu*((U.u(k,l)-U.u(k,l-1))/integ.dy+(interp_v_plus-interp_v_minus)/integ.dx);
    }
  }

  // * top left corner
  interp_rho = (U.rho(1,integ.ngy-1)+U.rho(1,integ.ngy-2))/2;
  mu = interp_rho*flow.nu;
  sigma(1,integ.ngy-1) = mu*(U.u(1,integ.ngy-1)-U.u(1,integ.ngy-2))/integ.dy;

  //* top right corner
  interp_rho = (U.rho(integ.ngx-2,integ.ngy-1)+U.rho(integ.ngx-2,integ.ngy-2))/2;
  mu = interp_rho*flow.nu;
  sigma(integ.ngx-2,integ.ngy-1) = mu*(U.u(integ.ngx-2,integ.ngy-1)-U.u(integ.ngx-2,integ.ngy-2))/integ.dy;


  return sigma;
}

Eigen::ArrayXXd sig_west(struct flowParams flow, struct integParams integ, struct flowQuant U) {

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  double interp_rho;
  double mu;

  double interp_u_plus;
  double interp_u_minus;

  int k,l;

  #pragma omp parallel for private(k,interp_rho,mu,interp_u_plus,interp_u_minus) collapse(2)
  for(int l = 1; l < integ.ngy-1; ++l) {
    for(int k = 1; k < integ.ngx; ++k) {
      interp_rho = (U.rho(k-1,l)+U.rho(k,l))/2;
      mu = interp_rho*flow.nu;
      interp_u_plus = (U.u(k,l)+U.u(k-1,l)+U.u(k-1,l+1)+U.u(k,l+1))/4;
      interp_u_minus = (U.u(k,l)+U.u(k,l-1)+U.u(k-1,l-1)+U.u(k-1,l))/4;

      sigma(k,l) = mu*((interp_u_plus-interp_u_minus)/integ.dy+(U.v(k,l)-U.v(k-1,l))/integ.dx);
    }
  }

  // ! handle top corners separately

  // * top left
  interp_rho = (U.rho(0,integ.ngy-2)+U.rho(1,integ.ngy-2))/2;
  mu = interp_rho*flow.nu;
  sigma(1,integ.ngy-2) = mu*(U.v(1,integ.ngy-2)-U.v(0,integ.ngy-2))/integ.dx;

  // * top right
  interp_rho = (U.rho(integ.ngx-1,integ.ngy-2)+U.rho(integ.ngx-2,integ.ngy-2))/2;
  mu = interp_rho*flow.nu;
  sigma(integ.ngx-1,integ.ngy-2) = (U.v(integ.ngx-1,integ.ngy-2)-U.v(integ.ngx-2,integ.ngy-2))/integ.dx;


  return sigma;
}

double temp(struct flowParams flow, double et, double u, double v) {

  // utility function for calculating the temperature at a point given a total
  // energy and velocity

  // et = 1/(gamma-1)*RT + 0.5*|V|^2
  // => (et-0.5*|V|^2)*(gamma-1)/R = T

  double T;

  T = (flow.gamma-1)/flow.R*(et-0.5*(u*u+v*v));

  return T;
}
