#pragma once

#ifndef VCP_LDBASE_FAST_FIXED_MINMAX_HPP
#define VCP_LDBASE_FAST_FIXED_MINMAX_HPP

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
typedef double AppData;
typedef vcp::pdblas AppPolicy;

typedef kv::interval< double > VData;
typedef vcp::pidblas VPOLICY;

typedef kv::dd ResData;
typedef vcp::mats< kv::dd > ResPOLICY;

typedef kv::interval< kv::dd > ResVData;
typedef vcp::imats< kv::dd > ResVPOLICY;

namespace vcp {
	template < typename _T, typename _P >
	struct VanDerPol : public vcp::Newton< _T, _P > {
		_T alpha;
		_T gamma;
		_T tau;
		_T beta;
		_T omega;
		_T mu;
		_T k;
		int order, vec_size;

		vcp::fourier_series< _T > x, y;
		vcp::fourier_series< _T > Bcos_nt;
		vcp::fourier_basis< _T, _P > Generator1, Generator2;

		void setting_newton( vcp::matrix< _T, _P >& zh ) override {
			x = Generator1.omit_vec_to_fourier_series( zh.submatrix( {0, vec_size - 1}, {0}) );
			y = Generator1.omit_vec_to_fourier_series( zh.submatrix( {vec_size, 2*vec_size - 1}, {0}) );

			Generator1.clear_data();
			Generator2.clear_data();
		}

		vcp::fourier_series< _T > func1(){
			return x.diff() - y/omega;
		}

		vcp::fourier_series< _T > func2(){
			return y.diff() - (-k*(x*x-1)*y - mu*x - gamma*x*x*x + alpha*x.delay(-omega*tau) + Bcos_nt)/omega;
		}

		vcp::matrix< _T, _P > f() override {
			Generator1.add_fourier_fx( this->func1() );
			Generator2.add_fourier_fx( this->func2() );
			vcp::matrix< _T, _P > zh = vercat( Generator1.output_fx(), Generator2.output_fx() );
			return zh;
		}

		vcp::matrix< _T, _P > Q1() {
			Generator1.clear_data();
			Generator2.clear_data();

			Generator1.add_dpt();

			Generator2.add_scalar_pt( -1/omega );

			return horzcat( Generator1.output_Jacobi(), Generator2.output_Jacobi() );
		}

		vcp::matrix< _T, _P > Q2() {
			Generator1.clear_data();
			Generator2.clear_data();

			Generator1.add_fourier_pt( 2*k*x*y/omega );
			Generator1.add_scalar_pt( mu/omega );
			Generator1.add_fourier_pt( 3/omega*gamma*x*x );
			Generator1.add_scalar_pt_delay( -alpha/omega, -omega*tau );

			Generator2.add_dpt();
			Generator2.add_fourier_pt( k*(x*x-1)/omega );

			return horzcat( Generator1.output_Jacobi(), Generator2.output_Jacobi() );
		}

    //ヤコビ行列を作る
		vcp::matrix< _T, _P > Df() override {
			return vercat( this->Q1(), this->Q2() );
		}

		template < typename _TM >
		void set_parameter( _TM Alpha, _TM Gamma, _TM Tau, _TM Beta, _TM Mu, _TM K, _TM Omega, int N){
			alpha = _T( Alpha );
			gamma = _T( Gamma );
			tau   = _T( Tau );
			beta  = _T( Beta );
			omega = _T( Omega )/_T(N);
			mu    = _T( Mu );
			k     = _T( K );

			Bcos_nt.zeros(N+1);
			Bcos_nt.set_cosm( _T( Beta ), N);
		}

		void set_order( const int& Order ){
			order = Order;
			this->x.zeros(order);
			this->y.zeros(order);
			this->Generator1.setting_m(order);
			this->Generator1.create_full_list();
			this->Generator2.setting_m(order);
			this->Generator2.create_full_list();
			vec_size = Generator1.list_size();
		}

		vcp::matrix< AppData, AppPolicy > FourierCoefficient(vcp::matrix< AppData, AppPolicy > u, vcp::matrix< AppData, AppPolicy > t, double tau, double T, int n,int m){

		  using std::cos;
		  using std::sin;

		  double h = 1./4096;
			int N = 10*1024;
		  int L = tau*4096;
			int l = 1000*N;
		  int ll = l/2;
			double pi = kv::constants<double>::pi();
		  double w1 = 2*pi*n/T;
		  int nn = w1/h+1;

		  vcp::matrix< AppData, AppPolicy > x, s;
		  x.zeros(nn+1, 1);
		  s.zeros(nn+1, 1);
			double a = 0;
		
		  for(int i=0; i<nn+1; i++){
		    x(i) = u(i+ll);
		    s(i) = 2.*pi*t(i+ll)/w1;
				a += x(i);
			}
		  a /= nn;

		  vcp::matrix< AppData, AppPolicy > f;
		  f.zeros(vec_size, 1);
		  f(0) = a;
		  double a1 = 0;
		  double b1 = 0;
		  for(int i=1; i<=(vec_size-1)/2; i++){
		    int j = 2*i;
		    a1 = 0;
		    b1 = 0;
		    for(int k=0; k<nn+1; k++){
		        a1 += x(k)*sin(i*s(k));
		        b1 += x(k)*cos(i*s(k));
		    }
		    a1 *= 2.;
		    a1 /= (double)nn;
		    b1 *= 2.;
		    b1 /= (double)nn;
		    f(j-1) = a1;
		    f(j) = b1;
		  }
		  return f;

		}

		double ff(double y,double u,double z,double t,double T){
			return -(k*(y*y-1)*u+mu*y+gamma*y*y*y-alpha*z-beta*cos(t*T));
		}
		vcp::matrix< AppData, AppPolicy > euler_vdp(double alpha, double gamma, double tau, double beta, double omega, int vec_size, int n){

		  vcp::matrix< AppData, AppPolicy > x;

		  std::ofstream result_file("euler.csv");

		  const double init_val = 0.0;


		  using std::cos;

			double T = omega;//T=1
		  double h = 1./4096;
		  int N = 10*1024;
		  int L = tau*4096;
			x.zeros(1000*N, 1);
		  double u = init_val;
		  double t,z,y;

		  for(int i=1; i<=(998*N-1); i++){
		      t=i*h;
		      if(i-L<=0){
						z = init_val;
					}else{
						z = x(i-1-L);
					}
					y = x(i-1);
					u = u + ff(y,u,z,t,T)*h;
					y = y+u*h;
					x(i) = y;
			}

			vcp::matrix< AppData, AppPolicy > t_;
		  t_.zeros(1000*N, 1);
			for(int i=0; i<1000*N; i++){
		      t_(i) = (i+1)*h;
		  }

			int l = 1000*N;
		  int ll = l/2;

			for(int i=1; i<=102000*5; i++){
				double uu = (x(ll+i)-x(ll+i-1))/h;
				result_file << x(ll+i-1) << "," << uu << std::endl;
			}

			return FourierCoefficient(x, t_, tau, T, n, vec_size);
		}

		vcp::matrix< AppData, AppPolicy > runge_kutta(double alpha, double gamma, double tau, double beta, double omega, int vec_size, int n){
			vcp::matrix < AppData, AppPolicy > x;
			std::ofstream runge("runge_kutta.csv");

			const double init_val = 0.0;

			using std::cos;

			double T = omega;//T=1
			double h = 1./4096;
			int N = 10*1024;
			int L = tau*4096;
			x.zeros(1000*N, 1);
			double y = init_val;
			double t,z1,z2,xx,yy;
			double kx1,kx2,kx3,kx4, ky1,ky2,ky3,ky4;
			for(int i=1; i<=(998*N-1); i++){
		        t=i*h;
		        if(i-L<=0){
					z1 = z2 = init_val;
				}else{
					z1 = x(i-1-L);
					z2 = x(i-L);
				}
				xx = x(i-1); //x(t)を保存
				yy = y; //x(t)を保存
// ルンゲクッタ法の公式
//      kx1 = h*fx(t,x,y)
//      kx2 = h*fx(t+h/2,x+kx1/2,y+ky1/2)
//      kx3 = h*fx(t+h/2,x+kx2/2,y+ky2/2)
//      kx4 = h*fx(t+h,x+kx3,y+ky3)

//      ky1 = h*fy(t,x,y)
//      ky2 = h*fy(t+h/2,x+kx1/2,y+ky1/2)
//      ky3 = h*fy(t+h/2,x+kx2/2,y+ky2/2)
//      ky4 = h*fy(t+h,x+kx3,y+ky3)
//       遅延項のところはn+(kx1/2)/h-Lなどを確かめないといけないので全てzで代用する.
//zをz+kx1/2などにする
				kx1 = yy * h;
				ky1 = ff(xx,yy,z1,t,T) * h;

				kx2 = ( yy + ky1/2 ) * h;
				ky2 = ff( xx + kx1/2, yy + ky1/2, (z1+z2)/2,t+h/2,T)*h;

				kx3 = ( yy + ky2/2 )*h;
				ky3 = ff( xx + kx2/2 , yy + ky2/2, (z1+z2)/2,t+h/2,T)*h;

				kx4 = ( yy + ky3/2 )*h;
				ky4 = ff( xx + kx3/2 , yy + ky3/2, z2 ,t,T)*h;

				x(i) = xx + ( kx1 + 2*kx2 + 2*kx3 + kx4 )/6;
				y = yy + ( ky1 + 2*ky2 + 2*ky3 + ky4 )/6;
			}
			vcp::matrix< AppData, AppPolicy > t_;
			t_.zeros(1000*N, 1);
			for(int i=0; i<1000*N; i++){
				t_(i) = (i+1)*h;
			}

			int l = 1000*N;
			int ll = l/2;

			for(int i=1; i<=102000*5; i++){
				double uu = (x(ll+i)-x(ll+i-1))/h;
				runge << x(ll+i-1) << "," << uu << std::endl;
			}

			return FourierCoefficient(x, t_, tau, T, n, vec_size);
		}
		vcp::matrix< AppData, AppPolicy > hein(double alpha, double gamma, double tau, double beta, double omega, int vec_size, int n){
			vcp::matrix < AppData, AppPolicy > x;
			std::ofstream hein("hein.csv");

			const double init_val = 0.0;

			using std::cos;

			double T = omega;//T=1
			double h = 1./4096;
			int N = 10*1024;
			int L = tau*4096;
			x.zeros(1000*N, 1);
			double y = init_val;
			double t,z1,z2,xx,yy;
			double kx1,kx2,kx3,kx4, ky1,ky2,ky3,ky4;
			for(int i=1; i<=(998*N-1); i++){
		        t=i*h;
		        if(i-L<=0){
					z1 = init_val;
					z2 = init_val;
				}else{
					z1 = x(i-1-L);
					z2 = x(i-1-L+1);
				}
				xx = x(i-1); //x(t)を保存
				yy = y; //x(t)を保存
// ホイン法の公式
//      kx1 = h*fx(t,x,y)
//      kx2 = h*fx(t+h/2,x+kx1/2,y+ky1/2)

//      ky1 = h*fy(t,x,y)
//      ky2 = h*fy(t+h/2,x+kx1/2,y+ky1/2)
//       遅延項のところはn+(kx1/2)/h-Lなどを確かめないといけないので全てzで代用する.
				kx1 = yy * h;
				ky1 = ff(xx,yy,z1,t,T) * h;

				kx2 = ( yy + ky1 ) * h;
				ky2 = ff( xx + kx1, yy + ky1, z2, t+h, T )*h;

				x(i) = xx + ( kx1 + kx2)/2;
				y = yy + ( ky1 + ky2 )/2;
			}
			vcp::matrix< AppData, AppPolicy > t_;
			t_.zeros(1000*N, 1);
			for(int i=0; i<1000*N; i++){
				t_(i) = (i+1)*h;
			}

			int l = 1000*N;
			int ll = l/2;

			for(int i=1; i<=102000*5; i++){
				double uu = (x(ll+i)-x(ll+i-1))/h;
				hein << x(ll+i-1) << "," << uu << std::endl;
			}
			return FourierCoefficient(x, t_, tau, T, n, vec_size);
		}

		vcp::matrix< AppData, AppPolicy > initial_fc(int n){
			vcp::matrix< AppData, AppPolicy > F;
		    //F = euler_vdp(alpha, gamma, tau, beta, omega*n, vec_size, n);
			F = runge_kutta(alpha, gamma, tau, beta, omega*n, vec_size, n);
			//F = hein(alpha, gamma, tau, beta, omega*n, vec_size, n);
			return F;
		}

		//void disp_continue() override {}
	};


	template < typename _T, typename _AccT, typename _P, typename _AccP >
	struct ResIter_VanDerPol : public vcp::Newton< _T, _P > {
		_T alpha;
		_T gamma;
		_T tau;
		_T beta;
		_T omega;
		_T mu;
		_T k;

		_AccT Acc_alpha;
		_AccT Acc_gamma;
		_AccT Acc_tau;
		_AccT Acc_beta;
		_AccT Acc_omega;
		_AccT Acc_mu;
		_AccT Acc_k;
		int order, vec_size;

		vcp::fourier_series< _T > x, y;
		vcp::fourier_series< _AccT > Acc_x, Acc_y;
		vcp::fourier_series< _AccT > Bcos_nt;
		vcp::fourier_basis< _T, _P > Generator1, Generator2;
		vcp::fourier_basis< _AccT, _AccP > AccGenerator1, AccGenerator2;


		void setting_newton( vcp::matrix< _T, _P >& zh ) override {
			x = Generator1.omit_vec_to_fourier_series( zh.submatrix( {0, vec_size - 1}, {0}) );
			y = Generator2.omit_vec_to_fourier_series( zh.submatrix( {vec_size, 2*vec_size - 1}, {0}) );

			vcp::matrix< _AccT > acc_zh;
			convert( zh, acc_zh );
			Acc_x = AccGenerator1.omit_vec_to_fourier_series( acc_zh.submatrix( {0, vec_size - 1}, {0}) );
			Acc_y = AccGenerator2.omit_vec_to_fourier_series( acc_zh.submatrix( {vec_size, 2*vec_size - 1}, {0}) );

			Generator1.clear_data();
			Generator2.clear_data();
			AccGenerator1.clear_data();
			AccGenerator2.clear_data();
		}

		vcp::fourier_series< _AccT > func1(){
			return Acc_x.diff() - Acc_y/Acc_omega;
		}

		vcp::fourier_series< _AccT > func2(){
			return Acc_y.diff() - ( -Acc_k*(Acc_x*Acc_x-1)*Acc_y - Acc_mu*Acc_x - Acc_gamma*Acc_x*Acc_x*Acc_x + Acc_alpha*Acc_x.delay(-Acc_omega*Acc_tau) + Bcos_nt)/Acc_omega;
		}

		vcp::matrix< _T, _P > f() override {
			AccGenerator1.add_fourier_fx( this->func1() );
			AccGenerator2.add_fourier_fx( this->func2() );
			vcp::matrix< _T, _P > zh;
			convert( vercat( AccGenerator1.output_fx(), AccGenerator2.output_fx() ), zh );
			return zh;
		}

		vcp::matrix< _T, _P > Q1() {
			Generator1.clear_data();
			Generator2.clear_data();

			Generator1.add_dpt();

			Generator2.add_scalar_pt( -1/omega );

			return horzcat( Generator1.output_Jacobi(), Generator2.output_Jacobi() );
		}

		vcp::matrix< _T, _P > Q2() {
			Generator1.clear_data();
			Generator2.clear_data();

			Generator1.add_dpt();
			Generator1.add_scalar_pt( k*(x*x-1)/omega );

			Generator2.add_scalar_pt( mu/omega );
			Generator2.add_fourier_pt( 3*gamma/omega*x*x );
			Generator2.add_scalar_pt_delay( alpha/omega, -omega*tau );

			return horzcat( Generator2.output_Jacobi(), Generator1.output_Jacobi() );
		}

		vcp::matrix< _T, _P > Df() override {
			return vercat( this->Q1(), this->Q2() );
		}

		template < typename _TM >
		void set_parameter( _TM Alpha, _TM Gamma, _TM Tau, _TM Beta, _TM Mu, _TM K, _TM Omega, int N){
			alpha = _T( Alpha );
			gamma = _T( Gamma );
			tau   = _T( Tau );
			mu    = _T( Mu );
			k     = _T( K );

			omega = _T(Omega)/_T(N);

			Acc_alpha = _AccT( Alpha );
			Acc_gamma = _AccT( Gamma );
			Acc_tau   = _AccT( Tau );
			Acc_mu    = _AccT( Mu );
			Acc_k     = _AccT( k );

			Acc_omega = _AccT(Omega)/_AccT(N);

			Bcos_nt.zeros(N+1);
			Bcos_nt.set_cosm( _AccT( Beta ), N);
		}

		void set_order( const int& Order ){
			order = Order;
			this->x.zeros(order);
			this->y.zeros(order);

			this->Generator1.setting_m(order);
			this->Generator1.create_full_list();

			this->Generator2.setting_m(order);
			this->Generator2.create_full_list();

			this->AccGenerator1.setting_m(order);
			this->AccGenerator1.create_full_list();

			this->AccGenerator2.setting_m(order);
			this->AccGenerator2.create_full_list();

			vec_size = Generator1.list_size();
		}

		//void disp_continue() override {}
	};

	template < typename _T, typename _P >
	_T bdab( const vcp::matrix< _T, _P >& G, const int m_app, const int m_verify, const int var_num = 1 ){
		int sm = (m_app*2 + 1)*var_num;
		// int sm = (m_app*4 + 1)*var_num;
		int mm = (m_verify*2 + 1)*var_num;
		vcp::matrix< _T, _P > A, B, CDf, D;
		A = G.submatrix( {0, sm-1},  {0, sm-1}  );
		B = G.submatrix( {0, sm-1},  {sm, mm-1} );
		CDf = G.submatrix( {sm, mm-1}, {0, sm-1}  );
		D = G.submatrix( {sm, mm-1}, {sm, mm-1} );

		vcp::matrix< _T, _P > TT, E, LL;
		TT = ltransmul(A);
	    LL.eye( TT.rowsize() );
	    LL(0,0) = 0.5;
		LL(1,1) = 0.5;
		eigsym(TT,E);
	    E = diag(E);
		// std::cout << "min(lambda) = " << min(E) << std::endl;
		std::cout << "\nCheck BDAB || DFinv ||_2 in " << VData(1)/sqrt((min(E))(0)) << std::endl;

		eigsymge(TT, LL, E);
		E = diag(E);
		std::cout << "\nCheck BDAB || DFinv L^1/2 ||_2 in " << VData(1)/sqrt((min(E))(0)) << std::endl;

		vcp::matrix< _T, _P > Dd, Df;
		Df = D;		
		Dd.zeros( mm-sm, mm-sm );
		Dd(0, 0) = D(0, 0);
		Dd(0, 1) = D(0, 1);

		Df(0, 0) = _T(0);
		Df(0, 1) = _T(0);
		for (int i = 1; i < mm-sm-1; i++){
			Dd(i, i-1) = D(i, i-1);
			Dd(i, i)   = D(i, i);
			Dd(i, i+1) = D(i, i+1);

			Df(i, i-1) = _T(0);
			Df(i, i) = _T(0);
			Df(i, i+1) = _T(0);
            
            //もしかしてこっちが正しい?
			//Dd(i, i) = D(i, i);
			//Dd(i, i+1) = D(i, i+1);
			//Dd(i+1, i) = D(i+1, i);
			//Dd(i+1, i+1)   = D(i+1, i+1);

			//Df(i, i) = _T(0);
			//Df(i, i+1) = _T(0);
			//Df(i+1, i) = _T(0);
			//Df(i+1, i+1) = _T(0);
		}
		Dd(mm-sm-1, mm-sm-2) = D(mm-sm-1, mm-sm-2);
		Dd(mm-sm-1, mm-sm-1) = D(mm-sm-1, mm-sm-1);
		
		Df(mm-sm-1, mm-sm-2) = _T(0);
		Df(mm-sm-1, mm-sm-1) = _T(0);
//		Df = D - Dd;

		vcp::matrix< _T, _P > invAB = abs(lss(A, B));
		_T kk1 = norminf( invAB )(0);
		std::cout << "(kk1): || A^{-1}*B ||_{inf} <= " << kk1 << std::endl;

		CDf = horzcat( CDf, Df );

		vcp::matrix< _T, _P > invDdCDf = lss( Dd, CDf );
		_T kk2 = norminf( invDdCDf )(0);
		std::cout << "(kk2): || Dd^{-1}*[C, Df] ||_{inf} <= " << kk2 << std::endl;

		_T kk = max(kk1,kk2);
		std::cout << "(kk): max( || A^{-1}*B ||_{inf}, || Dd^{-1}*[C, Df] ||_{inf}) <= " << kk << std::endl;

		_T norm_invA, norm_invDd;
		{
			vcp::matrix< _T, _P > invTMP, I;
			I.eye(A.rowsize());
			invTMP = lss(A, I);
			norm_invA = norminf(invTMP)(0);
			std::cout << "|| norm_invA || <= " << norm_invA << std::endl;

			I.clear();
			I.eye(D.rowsize());
			invTMP = lss(Dd, I);
			norm_invDd = norminf(invTMP)(0);
			std::cout << "|| Ddinv || <= " << norm_invDd << std::endl;
		}

		if ( kk.upper() >= 1 ){
			std::cout << "kk >= 1..." << std::endl;
			std::cout << "Verification Failure... (bijectivity is not verified: kk >= 1)" << std::endl;
			exit(0);
		}
		else{
			std::cout << "kk < 1 !!" << std::endl;
		}
		_T ninvM2 = norm_invA/(1 - kk);
		std::cout << "ninvM2 <= " << ninvM2 << std::endl;

		_T ninvM3 = norm_invDd/(1 - kk);
		std::cout << "ninvM3 <= " << ninvM3 << std::endl;

		_T Mn = max(ninvM2, ninvM3);
		return Mn;
	}

	std::vector< int > change_list( const int m_verify, const int var_num ){
		int mm1 = m_verify*2 + 1;
		int mm = mm1*var_num;

		std::vector< int > list;

		list.push_back( 0 );
		list.push_back( mm1 );
		for (int k = 1; k < mm1; k += 2){
			list.push_back( k );
			list.push_back( k+1 );
			list.push_back( mm1 + k );
			list.push_back( mm1 + k + 1 );
		}
		return list;
	}

	template < typename _T, typename _P >
	vcp::matrix< _T, _P > change_matrix( const vcp::matrix< _T, _P >& A, const std::vector< int >& clist ){
		vcp::matrix< _T, _P > B;
		B.zeros( A.rowsize(), A.columnsize() );

		std::cout << "A size: " << A.rowsize() << ", " << A.columnsize() << std::endl;
		std::cout << "clist size: " << clist.size() << std::endl;
		
		for(int i = 0; i < clist.size(); i++ ){
			int ii = clist[i];
			for(int j = 0; j < clist.size(); j++ ){
				int jj = clist[j];
				B(i, j) = A(ii, jj);
			}
		}
		
		return B;
	}


}

#endif
