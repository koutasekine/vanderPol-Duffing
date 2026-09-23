#include <iostream>
#include <fstream>
#include <cmath>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>


#include <vcp/fourier_basis.hpp>
#include <vcp/vcp_timer.hpp>

#include <vcp/newton.hpp>

#include "VanDerPol.hpp"

/* Approximate data type (Newton method) */
typedef double AppData;
typedef vcp::pdblas AppPolicy;

typedef kv::interval< double > VData;
typedef vcp::pidblas VPOLICY;

typedef kv::dd ResData;
typedef vcp::mats< ResData > ResPOLICY;

typedef kv::interval< ResData > ResVData;
typedef vcp::imats< ResData > ResVPOLICY;


int main(int argc, char *argv[]){

  std::cout.precision(17);
  if(argc!=10) {
    std::cout<<"./a.out [m_app] [alpha] [beta] [tomega] [tau] [gamma] [mu] [k] [n]"<<std::endl;
    exit(0);
  }

  int m_app = atoi(argv[1]); // 50
  int m_verify = 4*m_app;
  int app_list_size;
  vcp::matrix< AppData, AppPolicy > xh, yh, zh;
  vcp::matrix< VData, VPOLICY > ixh, iyh, izh;

  std::string alpha = argv[2]; // "1"
  std::string beta = argv[3]; // "5"
  std::string tomega = argv[4]; // "1"
  std::string tau = argv[5]; // "0.25"
  std::string gamma = argv[6]; // "1"
  std::string mu = argv[7]; // "0.1"
  std::string k = argv[8]; // "0.1"
  int n = atoi(argv[9]); // 1
  VData omega = VData(tomega)/VData(n);

  std::cout << "Create an approximate solution" << std::endl;
  {
    vcp::VanDerPol< AppData, AppPolicy > VDP;
    double pi = kv::constants< double >::pi();

    std::cout << "Approximate Fourier order m_app = " << m_app << std::endl;
    std::cout << "alpha = " << alpha << std::endl;
    std::cout << "beta = " << beta << std::endl;
    std::cout << "tilde{omega} = " << tomega << std::endl;
    std::cout << "omega = tilde{omega}/n = " << omega << std::endl;
    std::cout << "tau = " << tau << std::endl;
    std::cout << "gamma = " << gamma << std::endl;
    std::cout << "mu = " << mu << std::endl;
    std::cout << "k = " << k << std::endl;
    std::cout << "n = " << n << std::endl;

    VDP.set_order(m_app);
    VDP.set_parameter(
      std::stod(alpha),
      std::stod(gamma),
      std::stod(tau),
      std::stod(beta),
      std::stod(mu),
      std::stod(k),
      std::stod(tomega),
      n
    );
    // std::ofstream param("parameter.csv");
    // param << alpha << std::endl;
    // param << gamma << std::endl;
    // param << tau << std::endl;
    // param << beta << std::endl;
    // param << mu << std::endl;
    // param << k << std::endl;
    // param << omega << std::endl;
    // param << n << std::endl;
    // param.close();

    app_list_size = VDP.vec_size;
    xh.zeros(VDP.vec_size, 1);
    yh.zeros(VDP.vec_size, 1);

    std::string filebase="matrix_ds/before_xh_m"+std::to_string(m_app)+"_t"+tau+"_n"+std::to_string(n);
    // std::string filebase="before_xh_m"+std::to_string(m_app)+"_t"+tau+"_n"+std::to_string(n);
    std::fstream loadfile;
    loadfile.open(filebase + ".matrix_d", std::ios::in | std::ios::binary);
    if (loadfile.is_open()) {
      vcp::load(xh, filebase.c_str());
      std::cout<<"...loaded before xh data"<<std::endl;
    } else {

      xh = VDP.initial_fc(n); //こいつめちゃ時間かかる

    //  vcp::save(xh, filebase.c_str());
    //  std::cout<<"save xh data!"<<std::endl;
    }
    loadfile.close();

    // std::cout << "xh = \n" << xh << std::endl;

    for(int i=0; i<VDP.vec_size; i++){
      if(i == 0){
        yh(i) = 0;
      }else if(i%2 == 1){
        yh(i) = -xh(i+1)*(i/2+1);
      }else if(i%2 == 0){
        yh(i) = xh(i-1)*(i/2);
      }
    }

    // std::cout << "yh = \n" << yh << std::endl;

    zh = vercat( xh, yh );
    VDP.setting_newton( zh );
    std::ofstream beforensolv("beforensolv.csv");
    for(double w=0.; w<=n*pi; w+=0.001){
      beforensolv << VDP.x.value(w) << "," << (VDP.x.diff()/n).value(w) << "\n";
    }
    beforensolv.close();
    VDP.setting_newton_tol(32);
    zh = VDP.solve_nls( zh );
    VDP.setting_newton( zh );

    //start time
    vcp::time.tic();

    vcp::fourier_series< AppData > x = VDP.x;
    vcp::fourier_series< AppData > y = VDP.y;
    std::ofstream afternsolv("afternsolv.csv");
    for(double w=0.; w<=n*pi; w+=0.001){
      afternsolv << VDP.x.value(w) << "," << (VDP.x.diff()/n).value(w) << "\n";
    }
    afternsolv.close();
    //横軸t 縦軸x(t)
    std::ofstream xt("xt.csv");
    for(double w=0.; w<=n*pi; w+=0.001){
      xt << w*n << "," << VDP.x.value(w) << "\n";
    }
    xt.close();

    std::cout << "\n========================================================" << std::endl;
  	std::cout << "Approximation: " << std::endl;
  	// std::cout << "x = " << std::endl;
  	// std::cout << x << std::endl;

  	// std::cout << "y = " << std::endl;
  	// std::cout << y << std::endl;

    xh = zh.submatrix( {0, app_list_size - 1}, {0});
    yh = zh.submatrix( {app_list_size, 2*app_list_size - 1}, {0});
    {
      vcp::matrix< AppData, AppPolicy > tmp;
      tmp.zeros(m_verify*2+1,1);
      for(int i=0; i<xh.rowsize(); i++){
        tmp(i) = xh(i);
      }
      vcp::interval( tmp, ixh );
    }
    {
      vcp::matrix< AppData, AppPolicy > tmp;
      tmp.zeros(m_verify*2+1,1);
      for(int i=0; i<yh.rowsize(); i++){
        tmp(i) = yh(i);
      }
      vcp::interval( tmp, iyh );
    }
    izh = vercat(ixh, iyh);
  }

  std::cout << "\n========================================================" << std::endl;
  std::cout << "Compute the constant M" << std::endl;
  VData Mn, K0, M, xh_inf, yh_inf;
  {
    std::cout << "Verification Fourier order m_verify = " << m_verify << std::endl;
    VData pi = kv::constants< VData >::pi();
    vcp::VanDerPol< VData, VPOLICY > VDP;
    VDP.set_parameter(
      VData(alpha),
      VData(gamma),
      VData(tau),
      VData(beta),
      VData(mu),
      VData(k),
      VData(tomega),
      n
    );

    VDP.set_order(m_verify);
    VDP.setting_newton(izh);

    
    std::cout << "\nCompute || xh ||_{Linf} and || yh ||_{Linf}" << std::endl;
    vcp::fourier_series< VData > x = VDP.x;
    vcp::fourier_series< VData > y = VDP.y;
    xh_inf = x.Linfnorm(1000);
    yh_inf = y.Linfnorm(1000);
    std::cout << "|| xh ||_{Linf} = " << xh_inf << std::endl;
    std::cout << "|| yh ||_{Linf} = " << yh_inf << std::endl;

    { // K0 ( ||N'[zh]|| ) を計算
      vcp::matrix< VData, VPOLICY > K0_mat;
      K0_mat.zeros(2);
      std::cout << "\nCompute K0 ( ||N'[zh]|| )" << std::endl;
      K0_mat(0,0) = VData(0);
      K0_mat(0,1) = VData(1);
      K0_mat(1,0) = ( (2*VData(k)*(y) + 3*VData(gamma)*(x))*x +VData(mu)).Linfnorm(1000) + abs(VData(alpha));
      K0_mat(1,1) = VData(k)*(x*x-VData(1)).Linfnorm(1000);

      K0 = VData(1)/omega * normtwo(K0_mat)(0);

      std::cout << "K0 = " << K0 << std::endl;
    }


    vcp::matrix< VData, VPOLICY > DF = VDP.Df();
    VData sigman;
/*
    std::cout << "------------------------------ 1995's method -----------------------------" << std::endl;
    vcp::matrix< VData, VPOLICY > E, LL;
    DF = ltransmul(DF); //DF = (DF^t)*DF
    LL.eye( DF.rowsize() );
    LL(0,0) = 0.5;
    compsym(DF);
    eigsymge(DF, LL, E);

    E = diag(E);

    std::cout << "Eigenvalue problem: DF^T DF x = lambda x" << std::endl;
    std::cout << "All lambda" << std::endl;
    std::cout << "min(lambda) = " << min(E) << std::endl;
    std::cout << "max(lambda) = " << max(E) << std::endl;
    Mn = VData(1)/sqrt((min(E))(0));
    std::cout << "\nCheck Mn in " << Mn << std::endl;

    sigman = VData(1)/(VData(m_verify)+VData(1));
    std::cout << "sigman = " << sigman << std::endl;

    if(sigman*K0*(VData(1)+Mn*K0) >= 1){
      std::cout << "Error: verifyVDP: Compute M: sigman*K0(1+M_n*K0) >= 1..." << std::endl;
//      exit(0);
    }
    else{
      std::cout << "|| DF^{-1} ||_{L(Z,D(L))} <= (1+Mn*K0+sigman(1+Mn*K0)+Mn)/(1-sigman*K0(1+M_n*K0)) := M" << std::endl;
      M = (VData(1)+Mn*K0+sigman+(VData(1)+Mn*K0)+Mn)/(VData(1)-sigman*K0*(VData(1)+Mn*K0));
      std::cout << "M = " << M << std::endl;
    }

    DF = VDP.Df();      
    */
    // /*
    // std::ofstream testDF("testDF_94.csv");
    // for(int i=0;i<DF.rowsize();i++){
    //   for(int j=0;j<DF.rowsize();j++){
    //     if(j) testDF << ',' << DF(i,j);
    //     else testDF << DF(i,j);
    //   }
    //   testDF << std::endl;
    // }
    std::cout << "------------------------------ bdab method -----------------------------" << std::endl;
    std::vector< int > clist = vcp::change_list( m_verify, 2);
    vcp::matrix< VData, VPOLICY > cDF = change_matrix( DF, clist ); 

  //  std::cout << "DF = " << std::endl;
  //  std::cout << DF << std::endl;

  //  std::cout << "cDF = " << std::endl;
  //  std::cout << cDF << std::endl;


    // std::ofstream testcDF("testcDF_94.csv");
    // for(int i=0;i<cDF.rowsize();i++){
    //   for(int j=0;j<cDF.rowsize();j++){
    //     if(j) testcDF << ',' << cDF(i,j);
    //     else testcDF << cDF(i,j);
    //   }
    //   testDF << std::endl;
    // }
    // vcp::save(cDF,"testcDF_94"); 
    std::cout << "\n Mninf: " << std::endl;
		VData Mninf = bdab( cDF, m_app, m_verify, 2 );
		std::cout << "Mninf = " << Mninf << std::endl;

		std::cout << "\n Mnone: " << std::endl;
		VData Mnone = bdab( transpose(cDF), m_app, m_verify, 2 );
		std::cout << "Mnone = " << Mnone << std::endl;

		std::cout << "\n Mn: " << std::endl;
		Mn = sqrt(Mninf*Mnone);
		std::cout << "Mn = " << Mn << std::endl;

    sigman = VData(0);
    std::cout << "|| DF^{-1} ||_{L(Z,D(L))} <= (1+Mn*K0+sigman(1+Mn*K0)+Mn)/(1-sigman*K0(1+M_n*K0)) := M" << std::endl;
    //M = (VData(1)+Mn*K0+sigman*(VData(1)+Mn*K0)+Mn)/(VData(1)-sigman*K0*(VData(1)+Mn*K0));
    M = sqrt(pow(Mn,2) + pow(VData(1)+Mn*K0, 2));
    std::cout << "M = " << M << std::endl;

  }

  std::cout << "\n========================================================" << std::endl;
  std::cout << "Compute Residual norm || F(zh) ||_Z <= r " << std::endl;;
  VData r;
  {
    izh.clear();
    vcp::interval(zh, izh);
    vcp::VanDerPol< VData, VPOLICY > VDP;
    VDP.set_parameter(
      VData(alpha),
      VData(gamma),
      VData(tau),
      VData(beta),
      VData(mu),
      VData(k),
      VData(tomega),
      n
    );

    VDP.set_order(m_app);
    VDP.setting_newton(izh);

    vcp::fourier_series< VData > Fxh = VDP.func1();
    vcp::fourier_series< VData > Fyh = VDP.func2();
    VData r1 = Fxh.L2norm();
    VData r2 = Fyh.L2norm();
    
    r = sqrt( kv::constants< VData >::pi() ) * sqrt( pow(r1,2) + pow(r2,2) );

    r.lower() = r.upper();
    std::cout << "|| F(zh) ||_Z <= " << r << std::endl;

  }

  std::cout << "\n========================================================" << std::endl;
  VData eta, r0, br0;
  {
    std::cout << "Compute eta = M*r" << std::endl;
    eta = M*r;
    std::cout << "eta = " << eta << std::endl;

    std::cout << "Compute r0 = 2*eta = 2*M*r" << std::endl;
    r0 = VData(2)*eta;
    std::cout << "r0 = " << r0 << std::endl;

    std::cout << "Compute b(r0)" << std::endl;
    VData d = pow(VData(3), -0.25);
    VData tilde_br0, c1, c2;
    c1 = 2*VData(k)*(xh_inf+yh_inf+d*r0)+3*VData(gamma)*(2*xh_inf+d*r0);
    c2 = VData(k)*(2*xh_inf+d*r0);

    tilde_br0 = d*r0/omega * sqrt( pow(c1,2) + pow(c2,2) );

    std::cout << "tilde_b(r0) = " << tilde_br0 << std::endl;

    std::cout << "b(r0)" << std::endl;
    br0 = M*tilde_br0;
    std::cout << "b(r0) = " << br0 << std::endl;

    VData sufficient_cond_left = eta + r0*br0;
    std::cout << "eta + r0*br0 = " << sufficient_cond_left << std::endl;


    std::cout << "Sufficient condition" << std::endl;
    std::cout << sufficient_cond_left.upper() << " < " << r0.lower() << "?" << std::endl;

    //stop time
    vcp::time.toc();

    if (sufficient_cond_left.upper() < r0.lower()){
      std::cout << "Verification Success!!" << std::endl;
    }
    else{
      std::cout << "Verification Failure..." << std::endl;
    }

  }
}
