#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

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
#include <vcp/vcp_fio.hpp>

#include "improve_VanDerPol.hpp"

/* Approximate data type (Newton method) */
typedef double AppData;
typedef vcp::pdblas AppPolicy;

typedef kv::interval< double > VData;
typedef vcp::pidblas VPOLICY;

int main(int argc, char *argv[])
{
// ~~~~~~ パラメータを決定 
  std::cout.precision(17);
  if(argc!=10) {

    std::cout<<"./a.out [m_app] [alpha] [beta] [tomega] [tau] [gamma] [mu] [k] [n]"<<std::endl;

    exit(0);
  }

  int m_app = atoi(argv[1]);
  int m_verify = m_app*3 + 1;
  int app_list_size;
  vcp::matrix< AppData, AppPolicy > xh, yh, zh;
  vcp::matrix< VData, VPOLICY > ixhm, iyhm, izhm, ixh3m, iyh3m, izh3m;
  vcp::matrix< VData, VPOLICY > IDF_m;
  vcp::matrix< VData, VPOLICY > fm, DFm, Lm, Zm;
  vcp::matrix< VData, VPOLICY > N_3m_m, N_m_3m, Lbot, Cmat, Bmat;
  vcp::fourier_series< VData > xm;
  vcp::fourier_series< VData > ym;
  std::vector< int > clist = vcp::change_list(m_verify);
  std::vector< int > cmlist = vcp::change_list(m_app);
  VData xh_inf, yh_inf, sigmam = VData(1)/(VData(m_app)+VData(1)), K0, Mm, KmN, kappa, Sinv, CFinv, FinvB, Bnorm, etaZ, etaD, MZ,M;

  //文字列は文字データとして持つだけなので精度が悪くなることはない
  //別で変換するコードを設けることで誤差発生箇所を明示的にできる
  std::string alpha = argv[2];
  std::string beta = argv[3];
  std::string tomega = argv[4];
  std::string tau = argv[5];
  std::string gamma = argv[6];
  std::string mu = argv[7];
  std::string k = argv[8];
  int n = atoi(argv[9]);

  VData omega = VData(tomega)/VData(n);
// ~~~~~~~~~~~~~~~~~~


// ~~~~~~ 近似解を計算 ~~~~~
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
      std::stod(tomega), // set_parameter は内部で /n するので tilde{omega} を渡す
      n
    );
    app_list_size = VDP.vec_size;
    xh.zeros(VDP.vec_size, 1);
    yh.zeros(VDP.vec_size, 1);

    // ルンゲクッタから離散的な近似解を作る(別に精度はそれなりでOK)
    // 作った近似解から連続の関数を作る(フーリエ級数で表す)
    // 前にxhの計算後情報があればロードしてくる
    std::string filebase="matrix_ds/before_xh_m"+std::to_string(m_app)+"_t"+tau+"_n"+std::to_string(n);
    std::fstream loadfile;
    loadfile.open(filebase + ".matrix_d", std::ios::in | std::ios::binary);
    if (loadfile.is_open()) {
      vcp::load(xh, filebase.c_str());
      std::cout<<"...loaded before xh data"<<std::endl;
    } else {

      xh = VDP.initial_fc(n);

    //  vcp::save(xh, filebase.c_str());
    //  std::cout<<"save xh data!"<<std::endl;
    }
    loadfile.close();

    std::cout << "xh = \n" << xh << std::endl;

    // xからyの計算
    for(int i=0; i<VDP.vec_size; i++){
      if(i == 0){
        yh(i) = 0;
      }else if(i%2 == 1){
        yh(i) = -xh(i+1)*(i/2+1);
      }else if(i%2 == 0){
        yh(i) = xh(i-1)*(i/2);
      }
    }

    std::cout << "yh = \n" << yh << std::endl;

    zh = vercat( xh, yh );
    VDP.setting_newton( zh ); //newton法の初期値設定

    //ニュートン法計算前の結果出力
    std::ofstream beforensolv("beforensolv.csv");
    for(double w=0.; w<=n*pi; w+=0.001){
      beforensolv << VDP.x.value(w) << "," << (VDP.x.diff()/n).value(w) << "\n";
    }
    beforensolv.close();

    //得た関数は近似解を表すが,このままでは精度が悪い
    //このため,ガレルキン近似を利用してnewton法によりさらに精度のよい近似解を作成
    //newton法を利用して精度を上げる
    VDP.setting_newton_tol(1<<5);
    zh = VDP.solve_nls( zh );
    VDP.setting_newton( zh );

    std::cout<<"zh:"<<std::endl<<zh<<std::endl;

    //start time
    vcp::time.tic();

    //ここまでは係数だけの計算なのでフーリエ級数の形にする
    vcp::fourier_series< AppData > x = VDP.x;
    vcp::fourier_series< AppData > y = VDP.y;

    ////ニュートン法計算後の結果出力
    std::ofstream afternsolv("afternsolv.csv");
    for(double w=0.; w<=n*2*pi; w+=0.001){
      afternsolv << VDP.x.value(w) << "," << (VDP.x.diff()/n).value(w) << "\n";
    }
    afternsolv.close();

    std::cout << "\n========================================================" << std::endl;
  	std::cout << "Approximation: " << std::endl;
  	std::cout << "x = " << std::endl;
  	std::cout << x << std::endl;

  	std::cout << "y = " << std::endl;
  	std::cout << y << std::endl;

    vcp::interval(zh, izhm);

    std::cout << "Create matrix and vector" << std::endl;
    {
      vcp::VanDerPol< VData, VPOLICY > VDP;
      VDP.set_parameter(
        VData(alpha),
        VData(gamma),
        VData(tau),
        VData(beta),
        VData(mu),
        VData(k),
        VData(tomega), // set_parameter は内部で /n するので tilde{omega} を渡す
        n
      );

      VDP.set_order(m_app); // 精度保証のフーリエの次数
      VDP.setting_newton(izhm); // VDPにizhをもとにフーリエ級数の形を設定

      xm = VDP.x;
      ym = VDP.y;

      xh_inf = xm.Linfnorm(1000); //絶対値の最大を得る（引数は分割数）
      yh_inf = ym.Linfnorm(1000);

      fm = VDP.f();
      DFm = VDP.Df();
      Lm = VDP.make_L();
      Zm = lss(DFm,fm);

      {
        std::cout << "\nCompute Mm ( ||Fminv'|| )" << std::endl;        

        vcp::matrix<VData, VPOLICY> iLh; // L^{-1/2}（Lm は対角）
        iLh.zeros(Lm.rowsize(), Lm.columnsize());
        for(int i=0; i<Lm.rowsize(); i++) iLh(i,i) = 1/sqrt(Lm(i,i));
        auto G = ltransmul( iLh*DFm );                     // (L^{-1/2}F)^T (L^{-1/2}F) = F^T L^{-1} F

        vcp::matrix<VData, VPOLICY> lamda_E;
        eigsymge(G, Lm, lamda_E); //lamda_Lの固有値を抽出
        lamda_E = diag(lamda_E); //行列をベクトルにする

        Mm = VData(1)/sqrt((min(lamda_E))(0));
        std::cout << "Mm = " << Mm << std::endl;
      }

      // 区間型に変換
      ixhm = izhm.submatrix( {0, app_list_size - 1}, {0});
      iyhm = izhm.submatrix( {app_list_size, 2*app_list_size - 1}, {0});

      // m(m_app)からN(m_verify)に伸ばす
      {
        vcp::matrix< VData, VPOLICY > tmp;
        tmp.zeros(m_verify*2+1,1);
        for(int i=0; i < ixhm.rowsize(); i++){
          tmp(i) = ixhm(i);
        }
        ixh3m = tmp;
      }
      {
        vcp::matrix< VData, VPOLICY > tmp;
        tmp.zeros(m_verify*2+1,1);
        for(int i=0; i < iyhm.rowsize(); i++){
          tmp(i) = iyhm(i);
        }
        iyh3m = tmp;
      }
      izh3m = vercat(ixh3m, iyh3m);
    }
  }

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
      VData(tomega), // set_parameter は内部で /n するので tilde{omega} を渡す
      n
    );

    VDP.set_order(m_verify); //精度保証のフーリエの次数
    VDP.setting_newton(izh3m); //VDPにizhをもとにフーリエ級数の形を設定

    int up_m = (m_app*2+1)*2-1;
    int up_N = (m_verify*2+1)*2-1;
    
    {
      vcp::matrix<VData, VPOLICY> N = VDP.dN();
      N = vcp::change_matrix(N, clist);
      N_m_3m = N.submatrix({0, up_m}, {up_m+1, up_N});
      N_3m_m = N.submatrix({up_m+1, up_N}, {0, up_m});

      // L_botの計算
      vcp::matrix<VData, VPOLICY> L = VDP.L_bot();
      L = vcp::change_matrix(L, clist);
      Lbot = L.submatrix({up_m+1, up_N},{up_m+1, up_N});
      Bmat = N_m_3m;      
      Cmat = lss( Lbot, N_3m_m );
    }
  }

  vcp::matrix< VData, VPOLICY > DFs = vcp::change_matrix(DFm, cmlist);
  vcp::matrix< VData, VPOLICY > Ls  = vcp::change_matrix(Lm,  cmlist);
  vcp::matrix< VData, VPOLICY > Lhs; // Ls^{1/2}
  Lhs.zeros(Ls.rowsize(), Ls.columnsize());
  for(int i=0; i<Ls.rowsize(); i++) Lhs(i,i) = sqrt(Ls(i,i));


  { // K0 ( ||N'[zh]|| ) を計算
    std::cout << "\nCompute K0 ( ||N'[zh]|| )" << std::endl;    
    vcp::matrix< VData, VPOLICY > K0_mat;
    K0_mat.zeros(2);
    K0_mat(0,0) = VData(0);
    K0_mat(0,1) = VData(1);
    K0_mat(1,0) = ( (2*VData(k)*(ym) + 3*VData(gamma)*(xm))*xm +VData(mu)).Linfnorm(1000) + abs(VData(alpha));
    K0_mat(1,1) = VData(k)*(xm*xm-VData(1)).Linfnorm(1000);

    K0 = VData(1)/omega*normtwo(K0_mat)(0);

    std::cout << "K0 = " << K0 << std::endl;
  }
  {
    std::cout << "\nCompute KmN ( || C Finv B || )" << std::endl;
    auto G = Cmat*lss( DFs, Bmat );
    vcp::matrix<VData, VPOLICY> lamda_E;
    G = ltransmul(G);
    eigsym(G, lamda_E); //lamda_Lの固有値を抽出
    lamda_E = diag(lamda_E); //行列をベクトルにする

    KmN = sqrt(max(lamda_E)(0));
    std::cout << "KmN = " << KmN << std::endl;
  }
  {
    std::cout << "\nCompute kappa" << std::endl;
    kappa = sigmam*K0 + KmN;
    std::cout << "sigmam = " << sigmam << std::endl;
    std::cout << "sigmam*K0 = " << sigmam*K0 << std::endl;
    std::cout << "kappa = " << kappa << std::endl;
    std::cout << "(Old) kappa = " << sigmam*K0*(1+K0*Mm) << std::endl;

    Sinv = 1/(1-kappa);

    std::cout << "|| Sinv ||  = " << Sinv << std::endl;

  }

  {
    std::cout << "\nCompute CFinv ( || C Finv || )" << std::endl;
    vcp::matrix<VData, VPOLICY> I;
    I.eye(DFm.rowsize()); 

    auto G = Cmat*lss( DFs,  Ls);
    vcp::matrix<VData, VPOLICY> lamda_E;
    G = ltransmul(G);
    eigsymge(G, Ls, lamda_E); //lamda_Lの固有値を抽出

    lamda_E = diag(lamda_E); //行列をベクトルにする

    CFinv = sqrt(max(lamda_E)(0));
    std::cout << "|| C Finv || = " << CFinv << std::endl;
  }

  {
    std::cout << "\nCompute FinvB ( || Finv B || )" << std::endl; 
    auto G = Lhs*lss( DFs,  Bmat);
    vcp::matrix<VData, VPOLICY> lamda_E;
    G = ltransmul(G);
    eigsym(G, lamda_E); //lamda_Lの固有値を抽出
    lamda_E = diag(lamda_E); //行列をベクトルにする

    FinvB = sqrt(max(lamda_E)(0));
    std::cout << "|| Finv B || = " << FinvB << std::endl;
  }

  {
    std::cout << "\nCompute Bnorm ( || B || )" << std::endl; 
    //auto G = lss(Lm, Bmat);
    auto G = Bmat;
    vcp::matrix<VData, VPOLICY> lamda_E;
    G = ltransmul(G);
    eigsym(G, lamda_E); //lamda_Lの固有値を抽出
    lamda_E = diag(lamda_E); //行列をベクトルにする

    Bnorm = sqrt(max(lamda_E)(0));
    std::cout << "|| B || = " << Bnorm << std::endl;
  }

  //残差( || F(zh) ||_Z, || (I - Pm)F(zh) ||_Z, || Dinv F(zh) ||_Z, || C Dinv F(zh) ||_Z  )を計算 (etaの計算で使用)
  std::cout << "\n========================================================" << std::endl;
  std::cout << "Compute Residual norms || F(zh) ||_Z, || (I - Pm)F(zh) ||_Z, || Dinv F(zh) ||_Z, || C Dinv F(zh) ||_Z " << std::endl;;
  VData r_Fzh, r_IPmFzh, r_DinvFzh, r_CDinvFzh;
  {
    {
      vcp::VanDerPol< VData, VPOLICY > VDP;
      VDP.set_parameter(
        VData(alpha),
        VData(gamma),
        VData(tau),
        VData(beta),
        VData(mu),
        VData(k),
        VData(tomega), // set_parameter は内部で /n するので tilde{omega} を渡す
        n
      );

      VDP.set_order(m_app);
      VDP.setting_newton(izhm);

      vcp::fourier_series< VData > Fxh = VDP.func1();
      vcp::fourier_series< VData > Fyh = VDP.func2();

      VData r1 = Fxh.L2norm();
      VData r2 = Fyh.L2norm();

      // fourier_series::L2norm() は ||.||_{L^2(0,2pi)}/sqrt(pi) を返すため sqrt(pi) 倍して本文のノルムに戻す
      // d = 3^{-1/4} は正規化しない H^1 ノルムの定数なので、ノルムを本文の定義にそろえる
      r_Fzh = sqrt( kv::constants< VData >::pi() ) * sqrt( pow(r1,2) + pow(r2,2) );

      std::cout << "|| F(zh) ||_Z <= " << r_Fzh << std::endl;

      
      Fxh.set_a0(0);
      Fyh.set_a0(0);

      for(int i=1; i<=m_app; i++) {
        Fxh.set_sinm(0, i);
        Fxh.set_cosm(0, i);
        Fyh.set_sinm(0, i);
        Fyh.set_cosm(0, i);
      }

      r1 = Fxh.L2norm();
      r2 = Fyh.L2norm();

      r_IPmFzh = sqrt( kv::constants< VData >::pi() ) * sqrt( pow(r1,2) + pow(r2,2) );

      std::cout << "|| (I - Pm)F(zh) ||_Z <= " << r_IPmFzh << std::endl;


      VDP.setting_newton(Zm);
      Fxh = VDP.x;
      Fyh = VDP.y;
      r1 = Fxh.L2norm();
      r2 = Fyh.L2norm();
      r_DinvFzh = sqrt( kv::constants< VData >::pi() ) * sqrt( pow(r1,2) + pow(r2,2) );

      std::cout << "|| DFinv F(zh) ||_Z <= " << r_DinvFzh << std::endl;
    }
    {
      vcp::VanDerPol< VData, VPOLICY > VDP;
      VDP.set_parameter(
        VData(alpha),
        VData(gamma),
        VData(tau),
        VData(beta),
        VData(mu),
        VData(k),
        VData(tomega), // set_parameter は内部で /n するので tilde{omega} を渡す
        n
      );

      VDP.set_order(m_verify);
      auto Zmsort = Zm = change_vector(Zm, cmlist);
      
      //
      auto cc = lss(Lbot, N_3m_m*Zmsort);

      cc = change_3mbotTO3m(cc, m_app, m_verify);
      cc = inverse_change_vector(cc, clist);
      VDP.setting_newton(cc);
      vcp::fourier_series< VData > Fxh = VDP.x;
      vcp::fourier_series< VData > Fyh = VDP.y;

      VData r1 = Fxh.L2norm();
      VData r2 = Fyh.L2norm();

      r_CDinvFzh = sqrt( kv::constants< VData >::pi() ) * sqrt( pow(r1,2) + pow(r2,2) );

      std::cout << "|| C DFinv F(zh) ||_Z <= " << r_CDinvFzh << std::endl;
    }
  }

  {
    std::cout << "\n========================================================" << std::endl;
    std::cout << "Compute Newton term norms || F'[zh]^-1 F(zh) ||_Z" << std::endl;;

    etaZ = sqrt(pow(r_DinvFzh + Sinv*FinvB*( r_CDinvFzh + sigmam*r_IPmFzh ),2) + pow(Sinv*(r_CDinvFzh + sigmam*r_IPmFzh) ,2));
    std::cout << "|| F'[zh]^-1 F(zh) ||_Z <= " << etaZ << std::endl;
    etaD = sqrt( pow(etaZ,2) + pow( K0*etaZ + r_Fzh ,2) );
    std::cout << "|| F'[zh]^-1 F(zh) ||_D <= " << etaD << std::endl;
  }

  {
    vcp::matrix< VData, VPOLICY > Finv_mat;
    Finv_mat.zeros(2);


    // N_3m_m, N_m_3m, Lbot, Cmat, Bmat;
    // sigmam, K0, Mm, KmN, kappa, Sinv, CFinv, FinvB, Bnorm, etaZ, etaD
    Finv_mat(0,0) = abs( Mm + Sinv*FinvB*CFinv);
    Finv_mat(0,1) = abs( sigmam*Sinv*FinvB);
    Finv_mat(1,0) = abs( Sinv*CFinv );
    Finv_mat(1,1) = abs( sigmam*Sinv );
    MZ = normtwo(Finv_mat)(0);
    M = sqrt( pow(MZ,2) + pow(1+K0*MZ,2) );
    std::cout << "|| DFinv ||_{B(Z)} <= " << MZ << std::endl;
    std::cout << "|| DFinv ||_{B(Z, D)} <= " << M << std::endl;    

  }

  std::cout << "\n========================================================" << std::endl;
  // 結果の出力
  // etaからKを計算
  VData r0, br0;
  {
    std::cout << "Compute r0 = 2*eta" << std::endl;
    r0 = VData(2)*etaD;
    std::cout << "r0 = " << r0 << std::endl;

    std::cout << "Compute b(r0)" << std::endl;
    VData d = pow(VData(3), -0.25);
    VData tilde_br0, c1, c2;
    c1 = 2*VData(k)*(xh_inf+yh_inf+d*r0)+3*VData(gamma)*(2*xh_inf+d*r0);
    c2 = VData(k)*(2*xh_inf+d*r0);
    
    tilde_br0 = d*r0/omega * sqrt( pow(c1,2) + pow(c2,2) );

    std::cout << "tilde_b(r0) = " << tilde_br0 << std::endl;

    br0 = M*tilde_br0;
    std::cout << "b(r0) = " << br0 << std::endl;

    // stop time
    vcp:: time.toc();

    //b(r0)<1/2で十分条件
    if (br0.upper() < 0.5){
      std::cout << "Verification Success!!" << std::endl;
      std::cout << "\n|| u* - u^ || <= " << etaD/(1 - br0) << std::endl;
//=====================================================================
      //   真の解の最小周期が 2*pi（(1.2) では 1/n 分数調波）であることの確認
      //   (1.6) の外力 cos(n t) の周期は 2*pi/n なので、beta != 0 のとき解の周期は 2*pi/n の整数倍で、
      //   2*pi を割り切るもの（2*pi*k/n, k は n の約数）に限られる。
      //   最小周期が 2*pi 未満なら、n の素因数 q のいずれかについて h = 2*pi/q も周期になる。
      //   よって各素因数 q について「z* は h = 2*pi/q 周期でない」を示せばよい。
      //   z* が h 周期なら，平行移動の等長性と ||.||_Z <= ||.||_D より
      //     || zh(.+h) - zh ||_Z <= || zh(.+h) - z*(.+h) ||_Z + || z* - zh ||_Z <= 2*rho
      //   となる。一方，Parseval の等式より
      //     || zh(.+h) - zh ||_Z^2 = pi * sum_j 4 sin^2(j h/2) (as_j^2 + ac_j^2 + ys_j^2 + yc_j^2)
      //   （定数項は平行移動で不変）。これが 2*rho を真に上回れば h 周期は否定される。
      {
        VData rho = etaD/(1 - br0);
        VData pi = kv::constants< VData >::pi();
        std::cout << "\nCheck minimal period" << std::endl;
        if ( VData(beta).lower() <= 0 && VData(beta).upper() >= 0 ){
          std::cout << "beta may be 0: the minimal-period argument requires beta != 0" << std::endl;
        }
        else if ( n == 1 ){
          std::cout << "n = 1: nothing to check (period 2*pi is the forcing period)" << std::endl;
        }
        else{
          bool all_ok = true;
          int nn = n;
          for (int q = 2; q <= nn; q++){
            if (nn % q != 0) continue;
            while (nn % q == 0) nn /= q;      // q は n の素因数
            VData h = 2*pi/VData(q);
            VData s = VData(0);
            for (int j = 1; j <= m_app; j++){
              VData w = pow(sin(VData(j)*h/VData(2)), 2);
              s += VData(4)*w*( pow(xm.get_sin(j),2) + pow(xm.get_cos(j),2)
                               + pow(ym.get_sin(j),2) + pow(ym.get_cos(j),2) );
            }
            VData shift_norm = sqrt(pi*s);   // || zh(.+h) - zh ||_Z
            std::cout << "q = " << q << " : || zh(.+2pi/q) - zh ||_Z >= " << shift_norm.lower()
                      << " , 2*rho <= " << (VData(2)*rho).upper() << std::endl;
            if ( shift_norm.lower() > (VData(2)*rho).upper() ){
              std::cout << "  -> z* is not (2pi/" << q << ")-periodic" << std::endl;
            }
            else{
              std::cout << "  -> cannot exclude period 2pi/" << q << std::endl;
              all_ok = false;
            }
          }
          if (all_ok){
            std::cout << "Minimal period of z* is 2*pi (1/" << n << " subharmonic of (1.2))" << std::endl;
          }
          else{
            std::cout << "Minimal period is NOT verified" << std::endl;
          }
        }
      }
//=====================================================================
    }
    else{
      std::cout << "Verification Failure..." << std::endl;
    }
  }
  std::cout << "\n========================================================" << std::endl;
}
