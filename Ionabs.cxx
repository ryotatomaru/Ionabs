// XSPEC model subroutine for absorption from highly ionized gass   
//---
// number of model parameters: 
//      0: IonID:Atomic number + the number of electrons, e.g. H-like Fe is 2601, He-like Fe is 2602
//      1: Ion column [1.0e18 cm^{-2}]      
//      2: Kinematic temperature [keV] 3: Redshift (minus: blueshift)
#include <cmath>
#include <iostream> 
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <tuple>
#include <xsTypes.h>
using namespace std;
//#include <functionMap.h>
//#include <XSFunctions/Utilities/FunctionUtility.h>

// prototype from tableInterpolate.cxx
const double Pi = 3.141592654;
const double C =2.99792458e10 ; //[cm/s]
const double H_planck = 6.62607004e-27;// [erg s]
const double H_planck_keV = 4.13566766e-18;//[keV s]
const double Erg2keV=6.2415064e8;
const double keV2erg =1.60218e-9;
const double Amu = 931.494013*1.0e3; // [keV/c^2]
const double Coeff = 0.02654009; // pi*e**2/(m_e*c) [cm^2/s]
//vector<double> Trans_E, F_ij, A_vec;
//vector<double> Edge_E, Sigma_0, Slope, Cut_E;
//int Ion_id;
map<int, double> A_mass = 
{
    {6, 12.0096}, {7, 14.00643}, {8, 15.99903}, 
    {10, 20.1797}, {11, 22.990},  {12, 24.304},
    {13, 26.982},  {14, 28.084}, {15, 30.973}, 
    {16, 32.065},  {17, 35.453}, {18,39.948},
    {19, 39.0983}, {20, 40.078}, {21, 47.867}, 
    {22, 47.867},  {23,50.9415},{24, 51.9961}, 
	{25, 54.938049},{26, 55.845}, {27, 58.933200}, 
	{28,58.6934}
};
map<int, string> E_name= 
{
    {6, "C"}, {7, "N"}, {8, "O"}, 
    {10, "Ne"}, {11, "Na"},  {12, "Mg"},
    {13, "Al"},  {14, "Si"}, {15, "P"}, 
    {16, "S"},  {17, "Cl"}, {18,"Ar"},
    {19, "K"}, {20, "Ca"}, {21, "Sc"}, 
    {22, "Ti"}, {23, "V"},{24, "Cr"}, 
	{25, "Mn"},{26, "Fe"}, {27, "Co"}, 
	{28, "Ni"}
};
map<int, tuple<vector<double>, vector<double>, vector<double> > > Datamap_line;
map<int, tuple<vector<double>, vector<double>, vector<double>, vector<double> > > Datamap_bf;
void txt2datamap(map<int, tuple<vector<double>, vector<double>, vector<double> > >& linemap, 
		map<int, tuple<vector<double>, vector<double>, vector<double>, vector<double> > >& bfmap);
//map<int, vector<double> > Data_values;
double Voigt (double xx, double sigma, double lg, int r);
void calc_line_abs(double& n_ion, const tuple<vector<double>, vector<double>, vector<double> >& line_tuple, const double& b_v,  const RealArray& ene_with_z, RealArray& tau_g, RealArray& tau_mid);
void calc_bf_abs(double& n_ion, const tuple<vector<double>, vector<double>, vector<double>, vector<double>  >& bf_tuple, const RealArray& ene_with_z, RealArray& tau_g, RealArray& tau_mid);
extern "C" {
    //void Kabs(const RealArray& energy, const RealArray& parameter,
            //int spectrum, RealArray& flux, RealArray& fluxError, 
            //const string& init)
	void Ionabs(const RealArray& energy, const RealArray& parameter, 
				   int spectrum, RealArray& flux, RealArray& fluxVariance,
				   const string& init){
		//energy [keV], flux  [photons/cm^2/s] 
		//RealArray eparams(4);
		int antr = (int) parameter[0];//Atomic number + number of electrons + quantum principal numter-1 (example, Fe XXV Kalpha= 2620)
		double N_z = parameter[1]*1.0e18;   //  Ion column [cm^{-2}]
		double kT = parameter[2]; //kT [keV]
		double z = parameter[3];// z 
		double doppler = 1.0+z;
		int Nflux=energy.size();
		//int num_a = 21;
		//double N_z_a[num_a*2];
		RealArray tau_grid(Nflux);
		RealArray tau_mid(Nflux-1);
		RealArray e_with_z(Nflux);
		flux.resize(Nflux-1);
		fluxVariance.resize(0);
		for (int j=0;j<Nflux;j++){
			e_with_z[j]=energy[j]*doppler;
			tau_grid[j]=0.0;
		}
		for (int j=0;j<Nflux-1;j++){
			tau_mid[j]=0.0;
			flux[j]=1.0;
		}
		if (Datamap_line.empty()){
			cout<< "Loading data files" <<endl;
			txt2datamap(Datamap_line, Datamap_bf);
			if (Datamap_line.size()==0){
				cout << "Data sets haven't load." <<endl;
				cout << "Something is wrong!" <<endl;
				exit(0);
			}
		}
		if (Datamap_line.count(antr) ==0) {
			cout <<"No data are available for this ion." <<endl;
			cout <<"The model calculates absorption by H-like Fe instead." <<endl;
			antr=2601;
		}
		/*
		vector<double> v1, v2, v3;
		v1= get<0>(Datamap_line[antr]);
		v2= get<1>(Datamap_line[antr]);
		v3= get<2>(Datamap_line[antr]);
		for (int i=0;i<v1.size(); i++){
			cout <<"E_t ="<< v1[i]<< ", f_ij =" <<v2[i] << ", Aij=" <<v3[i]<<endl;
		}

		for (auto it = Datamap_line.begin(); it != Datamap_line.end(); ++it) {
        std::cout << "Key: " << it->first << std::endl;
		}
		*/


		///loading atomic data sets///
		int a_num = antr/100;
		int num_e = antr%100;
	

		//cout << a_num<< endl;
		double atom_q = A_mass[a_num];
		double m_atom = atom_q*Amu;// keV/c^2
		//int num_id = A[antr];
		double b_v;
		if (kT>0.0) {
			b_v = sqrt(2.0*kT/m_atom); //'b' for the Voigt function  ( Dimensionless vth/c) 
			  //This returns effective temperature  for the Doppler broadening.
		}else{
			b_v = abs(kT); //'b' for the Voigt function  ( Dimensionless vth/c) 
		}
		/// end loading atomic data to ///
		calc_line_abs(N_z, Datamap_line[antr], b_v,  e_with_z, tau_grid, tau_mid);
		calc_bf_abs(N_z, Datamap_bf[antr], e_with_z, tau_grid, tau_mid );
		for (int i =0; i<Nflux-1; i++){
			if (tau_grid[i]<1.0e-5 && tau_mid[i]<1.0e-5) {
				flux[i]= 1.0;
			}else{
				flux[i]= 1.0/6.0*(exp(-tau_grid[i])+exp(-tau_grid[i+1])+4.0*exp(-tau_mid[i]));
			}
		}
		//double b_v = sqrt(2.0*kT/m_atom); //'b' for the Voigt function  ( Dimensionless vth/c) 
		return ;
	}
};
void txt2datamap(map<int, tuple<vector<double>, vector<double>, vector<double> > >& linemap, 
		map<int, tuple<vector<double>, vector<double>, vector<double>, vector<double> > >& bfmap){
//map<int, tuple<vector<double>, vector<double>, vector<double> > > Datamap_line;
	//string str_buf, str_conma_buf;
	//vector<double> vals;
	vector<double> trans_e, f_ij, a_ij, edge_e, sigma_0, slope, cut_e;
	const char* env_p = getenv("IONABS_DATA_PATH");
	if (env_p==nullptr){
		cout<< "The environment variable IONABS_DATA_PATH was not found" <<endl;
		cout<< "Please set that to data_ionabs" <<endl;
		exit(0);
	}
	string env_str = string(env_p);
	string element, num_e, f_name_line, f_name_bf;
	double num0, num1, num2, num3, num4, num5, num6, num7, num8;
	int key;
	for (auto iter=E_name.begin();iter !=E_name.end(); iter++) {
		element= iter->second;
		key=iter->first;
		//cout <<"Element is "<< element <<endl;
		for (int k=0; k<3;k++){
			if (k>8) {
				num_e =to_string(k+1);
			}else{
				num_e = "0"+to_string(k+1);
			}
			f_name_line= env_str+"/"+element+num_e+".px.tr";
			f_name_bf= env_str+"/"+element+num_e+".pi";
			//cout << "Line file name is " << f_name_line<<endl;
			//cout << "Edge file name is " << f_name_bf<<endl;
			ifstream ifs(f_name_line);
			cout << "Reading " << f_name_line<<endl;
			if (ifs){
				while(ifs >> num0 >> num1 >> num2 >> num3 >> num4 >> num5 >> num6 ){
					trans_e.push_back(num4*1.0e-3);
					f_ij.push_back(num5);
					a_ij.push_back(num6);
				}
			}else{
				cout << "Ion line data can not be found!" <<endl;
			}
			cout << "Finish reading" << f_name_line <<endl;
			//cout << "Key is " << key+k+1<<endl;
			linemap[100*key+k+1] = make_tuple(trans_e, f_ij, a_ij);
			trans_e.clear();
			f_ij.clear();
			a_ij.clear();
			

			ifstream ifs_b(f_name_bf);
			cout << "Reading " << f_name_bf<<endl;
			if (ifs_b){
				while(ifs_b >> num0 >> num1 >> num2 >> num3 >> num4 >> num5 >> num6 >> num7 >>num8){
					edge_e.push_back(num5*1.0e-3);
					sigma_0.push_back(num6);
					slope.push_back(num7);
					cut_e.push_back(num8*1.0e-3);
				}
			}else{
				cout << "Edge data can not be found!" <<endl;
			}
			cout << "Finish reading" << f_name_bf<<endl;
			//cout << "Key is " << key+k +1<<endl;
			bfmap[100*key+k+1] = make_tuple(edge_e, sigma_0, slope, cut_e);
			edge_e.clear();
			sigma_0.clear();
			slope.clear();
			cut_e.clear();
		}
    }
	return;
};
    

void calc_line_abs(double& n_ion, const tuple<vector<double>, vector<double>, vector<double> >& line_tuple, const double& b_v,  const RealArray& ene_with_z, RealArray& tau_g, RealArray& tau_mid){
	vector<double> t_ene, f_vec, A_vec;
	t_ene = get<0>(line_tuple);
	f_vec= get<1>(line_tuple);
	A_vec= get<2>(line_tuple);
	int tran_n= A_vec.size();
	int Nflux = ene_with_z.size();
	//vector<double> tau_g(Nflux,0.0); 
	//vector<double> tau_mid(Nflux-1,0.0); 
	double a_v, u_v, H_v, phi_v, tau_v, ene ;
	//int kk=0;
	//cout<< "The number of lines are " << tran_n<< endl;
	//tran_n = 2;
    for (int jj=0; jj<tran_n;jj++){
   //for (int jj=0; jj<1;jj++){
		a_v = A_vec[jj]*H_planck_keV/(4.0*Pi*t_ene[jj]*b_v);// 
		//cout << "a_v of ionabs is "<< a_v <<endl;
		//cout << "b_v of ionabs is "<< b_v <<endl;
        for (int i=0; i<Nflux; i++){
            ene = ene_with_z[i];
			u_v=0.0;
			H_v=0.0;
			tau_v = 0.0;
			phi_v = 0.0;
            if (ene > t_ene[jj]*0.5 and ene < t_ene[jj]*1.5){
				u_v = (ene - t_ene[jj])/(t_ene[jj]*b_v);
				H_v = sqrt(Pi)*Voigt(u_v, 1.0/sqrt(2.0), a_v, 5);
				//phi_v[jj] = 1.0/(E_c[jj]*b_v*sqrt(Pi))*H_v[jj];
				phi_v = H_planck_keV/(t_ene[jj]*b_v*sqrt(Pi))*H_v; //[s]
				//tau_v[jj] = 3.496e-20*N_z*Pi*ful[jj]*phi_v[jj];
				tau_v = n_ion*f_vec[jj]*Coeff*phi_v; //[cm^-2 cm^2/s s]
				tau_g[i]+=tau_v;
				//flux_g[i] =flux_g[i]*exp(-tau_v);
				//flux_g[i] *= exp(-tau_v);
				//cout << "u_v = "<< u_v <<endl;
			}else{
            tau_g[i] += 0.0;
			}
			
	   }
   }

	for (int jj=0; jj<tran_n;jj++){
		a_v = A_vec[jj]*H_planck_keV/(4.0*Pi*t_ene[jj]*b_v);// 
   		for (int i=0; i<Nflux-1; i++){
            //ene = energy[i]*(1.0+z);
            ene= (ene_with_z[i]+ene_with_z[i+1])*0.5;
			u_v=0.0;
			H_v=0.0;
			tau_v = 0.0;
			phi_v = 0.0;
            if (ene > t_ene[jj]*0.5 and ene < t_ene[jj]*1.5){
				u_v = (ene - t_ene[jj])/(t_ene[jj]*b_v);
				H_v = sqrt(Pi)*Voigt(u_v, 1.0/sqrt(2.0), a_v, 5);
				//phi_v[jj] = 1.0/(E_c[jj]*b_v*sqrt(Pi))*H_v[jj];
				phi_v = H_planck_keV/(t_ene[jj]*b_v*sqrt(Pi))*H_v; //[s]
				//tau_v[jj] = 3.496e-20*N_z*Pi*ful[jj]*phi_v[jj];
				tau_v = n_ion*f_vec[jj]*Coeff*phi_v; //[cm^-2 cm^2/s s]
				//cout << "tau_" << jj << tau_v[jj] <<endl;
				//flux_mid[i] = flux_mid[i]*exp(-tau_v);
				tau_mid[i] +=tau_v;
			}else{
             tau_mid[i] +=0.0;
			}
			
		}
   	}

	/*
	for (int i=0; i<Nflux;i++){
		
		if (tau_g[i] >1.0e-10){
			flux_g[i] = exp(-tau_g[i]);
		}else{
			flux_g[i] = 1.0;
		}
	   //cout << "tau_g of " << i << " is " <<tau_g[i] <<endl;
	}
	for (int i=0; i<Nflux-1;i++){
		if (tau_mid[i] >1.0e-10){
			flux_mid[i] = exp(-tau_mid[i]);
		}else{
			flux_mid[i] = 1.0;
		}
	   //cout << "tau_mid of " << i << " is " <<tau_mid[i] <<endl;
	}

	for (int i=0; i<Nflux-1;i++){
	   flux[i] *= 1.0/6.0*(flux_g[i]+4.0*flux_mid[i]+flux_g[i+1]);// numerical integration divided by bin size (Simpson's rule, flux/bin)
	   //flux[i] =(flux_g[i]+flux_g[i+1])*0.5;
    }
	*/
	return ;
};
void calc_bf_abs(double& n_ion, const tuple<vector<double>, vector<double>, vector<double>, vector<double> >& bf_tuple, const RealArray& ene_with_z, RealArray& tau_g, RealArray& tau_mid){
	//vector<double> edge_ene , sigma_0, slope, cut_e;
    vector<double> edge_e, sigma_0, slope, cut_e;
	edge_e = get<0>(bf_tuple);
	sigma_0 = get<1>(bf_tuple);
	slope = get<2>(bf_tuple);
	cut_e = get<3>(bf_tuple);
	int tran_n= edge_e.size();
	int Nflux = ene_with_z.size();
	double tau, ene;
	for (int j=0; j< tran_n;j++){
		for (int i=0; i<Nflux;i++){
           ene = ene_with_z[i];
		   if (ene>=edge_e[j]){
			tau= n_ion*sigma_0[j]*pow((ene/edge_e[j]), slope[j])*exp(-(ene/cut_e[j]));
		   } 
		   else{
			   tau=0.0;
		   }
		   tau_g[i] += tau;
		   //flux_g[i] = 1.0/6.0*(flux_grid[i]+4.0*flux_mid[i]+flux_grid[i+1]);// numerical integration divided by bin size (Simpson's rule, flux/bin)
		}

        for (int i=0; i<Nflux-1; i++){
            ene= (ene_with_z[i]+ene_with_z[i+1])*0.5;
		   if (ene>=edge_e[j]){
			tau= n_ion*sigma_0[j]*pow((ene/edge_e[j]), slope[j])*exp(-(ene/cut_e[j]));
			//tau= n_ion*sigma_0*pow((ene/edge_e), slope)*exp(-(ene/cue_e));
		   } 
		   else{
			   tau=0.0;
		   }
		   tau_mid[i] += tau;
		}
	}

	/*
	for (int i=0; i<Nflux-1;i++){
	   //delta_e =ene_with_z[i+1]-ene_with_z[i]
	   flux[i] *= 1.0/6.0*(flux_g[i]+4.0*flux_mid[i]+flux_g[i+1]);// numerical integration divided by bin size (Simpson's rule, flux/bin)
    }
	*/
	return ;
};


////////////////////////////////////////////////////////////////////////////////
/// Computation of Voigt function (normalised).
/// Voigt is a convolution of the two functions:
/// \f[
/// gauss(xx) = \frac{1}{(\sqrt{2\pi} sigma)} e^{\frac{xx^{2}}{(2 sigma{^2})}}
/// \f]
/// and
/// \f[
/// lorentz(xx) = \frac{ \frac{1}{\pi} \frac{lg}{2} }{ (xx^{2} + \frac{lg^{2}}{4}) }
/// \f]
/// .
///
/// The Voigt function is known to be the real part of Faddeeva function also
/// called complex error function [2].
///
/// The algoritm was developed by J. Humlicek [1].
/// This code is based on fortran code presented by R. J. Wells [2].
/// Translated and adapted by Miha D. Puc
///
/// To calculate the Faddeeva function with relative error less than 10^(-r).
/// r can be set by the the user subject to the constraints 2 <= r <= 5.
///
///  - [1] J. Humlicek, JQSRT, 21, 437 (1982).
///  - [2] [R.J. Wells "Rapid Approximation to the Voigt/Faddeeva Function and its Derivatives" JQSRT 62 (1999), pp 29-48.](http://www-atm.physics.ox.ac.uk/user/wells/voigt.html)
double Voigt (double xx, double sigma, double lg, int r)
{
   if ((sigma < 0 || lg < 0) || (sigma==0 && lg==0)) {
      return 0;  // Not meant to be for those who want to be thinner than 0
   }

   if (sigma == 0) {
      return lg * 0.159154943  / (xx*xx + lg*lg /4); //pure Lorentz
   }

   if (lg == 0) {   //pure gauss
      return 0.39894228 / sigma *exp(-xx*xx / (2*sigma*sigma));
   }

   double x, y, k;
   x = xx / sigma / 1.41421356;
   y = lg / 2 / sigma / 1.41421356;

   double r0, r1;

   if (r < 2) r = 2;
   if (r > 5) r = 5;

   r0=1.51 * exp(1.144 * (double)r);
   r1=1.60 * exp(0.554 * (double)r);

   // Constants

   const double rrtpi = 0.56418958;  // 1/SQRT(pi)

   double y0, y0py0, y0q;                      // for CPF12 algorithm
   y0 = 1.5;
   y0py0 = y0 + y0;
   y0q = y0 * y0;

   double c[6] = { 1.0117281, -0.75197147, 0.012557727, 0.010022008, -0.00024206814, 0.00000050084806};
   double s[6] = { 1.393237, 0.23115241, -0.15535147, 0.0062183662, 0.000091908299, -0.00000062752596};
   double t[6] = { 0.31424038, 0.94778839, 1.5976826, 2.2795071, 3.0206370, 3.8897249};

   // Local variables

   int j;                                        // Loop variables
   int rg1, rg2, rg3;                            // y polynomial flags
   double abx, xq, yq, yrrtpi;                 // --x--, x^2, y^2, y/SQRT(pi)
   double xlim0, xlim1, xlim2, xlim3, xlim4;   // --x-- on region boundaries
   double a0=0, d0=0, d2=0, e0=0, e2=0, e4=0, h0=0, h2=0, h4=0, h6=0;// W4 temporary variables
   double p0=0, p2=0, p4=0, p6=0, p8=0, z0=0, z2=0, z4=0, z6=0, z8=0;
   double xp[6], xm[6], yp[6], ym[6];          // CPF12 temporary values
   double mq[6], pq[6], mf[6], pf[6];
   double d, yf, ypy0, ypy0q;

   //***** Start of executable code *****************************************

   rg1 = 1;  // Set flags
   rg2 = 1;
   rg3 = 1;
   yq = y * y;  // y^2
   yrrtpi = y * rrtpi;  // y/SQRT(pi)

   // Region boundaries when both k and L are required or when R<>4

   xlim0 = r0 - y;
   xlim1 = r1 - y;
   xlim3 = 3.097 * y - 0.45;
   xlim2 = 6.8 - y;
   xlim4 = 18.1 * y + 1.65;
   if ( y <= 1e-6 ) {                      // When y<10^-6 avoid W4 algorithm
      xlim1 = xlim0;
      xlim2 = xlim0;
   }

   abx = fabs(x);                                // |x|
   xq = abx * abx;                               // x^2
   if ( abx > xlim0 ) {                          // Region 0 algorithm
      k = yrrtpi / (xq + yq);
   } else if ( abx > xlim1 ) {                   // Humlicek W4 Region 1
      if ( rg1 != 0 ) {                          // First point in Region 1
         rg1 = 0;
         a0 = yq + 0.5;                          // Region 1 y-dependents
         d0 = a0*a0;
         d2 = yq + yq - 1.0;
      }
      d = rrtpi / (d0 + xq*(d2 + xq));
      k = d * y * (a0 + xq);
   } else if ( abx > xlim2 ) {                   // Humlicek W4 Region 2
      if ( rg2 != 0 ) {                          // First point in Region 2
         rg2 = 0;
         h0 = 0.5625 + yq * (4.5 + yq * (10.5 + yq * (6.0 + yq)));
                                                 // Region 2 y-dependents
         h2 = -4.5 + yq * (9.0 + yq * ( 6.0 + yq * 4.0));
         h4 = 10.5 - yq * (6.0 - yq * 6.0);
         h6 = -6.0 + yq * 4.0;
         e0 = 1.875 + yq * (8.25 + yq * (5.5 + yq));
         e2 = 5.25 + yq * (1.0 + yq * 3.0);
         e4 = 0.75 * h6;
      }
      d = rrtpi / (h0 + xq * (h2 + xq * (h4 + xq * (h6 + xq))));
      k = d * y * (e0 + xq * (e2 + xq * (e4 + xq)));
   } else if ( abx < xlim3 ) {                   // Humlicek W4 Region 3
      if ( rg3 != 0 ) {                          // First point in Region 3
         rg3 = 0;
         z0 = 272.1014 + y * (1280.829 + y *
                              (2802.870 + y *
                               (3764.966 + y *
                                (3447.629 + y *
                                 (2256.981 + y *
                                  (1074.409 + y *
                                   (369.1989  + y *
                                    (88.26741 + y *
                                     (13.39880 + y)
                                     ))))))));   // Region 3 y-dependents
         z2 = 211.678 + y * (902.3066 + y *
                             (1758.336 + y *
                              (2037.310 + y *
                               (1549.675 + y *
                                (793.4273 + y *
                                 (266.2987 + y *
                                  (53.59518 + y * 5.0)
                                  ))))));
         z4 = 78.86585 + y * (308.1852 + y *
                              (497.3014 + y *
                               (479.2576 + y *
                                (269.2916 + y *
                                 (80.39278 + y * 10.0)
                                 ))));
         z6 = 22.03523 + y * (55.02933 + y *
                              (92.75679 + y *
                               (53.59518 + y * 10.0)
                               ));
         z8 = 1.496460 + y * (13.39880 + y * 5.0);
         p0 = 153.5168 + y * (549.3954 + y *
                              (919.4955 + y *
                               (946.8970 + y *
                                (662.8097 + y *
                                 (328.2151 + y *
                                  (115.3772 + y *
                                   (27.93941 + y *
                                    (4.264678 + y * 0.3183291)
                                    )))))));
         p2 = -34.16955 + y * (-1.322256+ y *
                               (124.5975 + y *
                                (189.7730 + y *
                                 (139.4665 + y *
                                  (56.81652 + y *
                                   (12.79458 + y * 1.2733163)
                                   )))));
         p4 = 2.584042 + y * (10.46332 + y *
                              (24.01655 + y *
                               (29.81482 + y *
                                (12.79568 + y * 1.9099744)
                                )));
         p6 = -0.07272979 + y * (0.9377051 + y *
                                 (4.266322 + y * 1.273316));
         p8 = 0.0005480304 + y * 0.3183291;
      }
      d = 1.7724538 / (z0 + xq * (z2 + xq * (z4 + xq * (z6 + xq * (z8 + xq)))));
      k = d * (p0 + xq * (p2 + xq * (p4 + xq * (p6 + xq * p8))));
   } else {                             // Humlicek CPF12 algorithm
      ypy0 = y + y0;
      ypy0q = ypy0 * ypy0;
      k = 0.0;
      for (j = 0; j <= 5; j++) {
         d = x - t[j];
         mq[j] = d * d;
         mf[j] = 1.0 / (mq[j] + ypy0q);
         xm[j] = mf[j] * d;
         ym[j] = mf[j] * ypy0;
         d = x + t[j];
         pq[j] = d * d;
         pf[j] = 1.0 / (pq[j] + ypy0q);
         xp[j] = pf[j] * d;
         yp[j] = pf[j] * ypy0;
      }
      if ( abx <= xlim4 ) {                      // Humlicek CPF12 Region I
         for (j = 0; j <= 5; j++) {
            k = k + c[j]*(ym[j]+yp[j]) - s[j]*(xm[j]-xp[j]) ;
         }
      } else {                                   // Humlicek CPF12 Region II
         yf = y + y0py0;
         for ( j = 0; j <= 5; j++) {
            k = k + (c[j] *
                 (mq[j] * mf[j] - y0 * ym[j])
                    + s[j] * yf * xm[j]) / (mq[j]+y0q)
                 + (c[j] * (pq[j] * pf[j] - y0 * yp[j])
                   - s[j] * yf * xp[j]) / (pq[j]+y0q);
         }
         k = y * k + exp(-xq);
      }
   }
   return k / 2.506628 / sigma; // Normalize by dividing by sqrt(2*pi)*sigma.
};

