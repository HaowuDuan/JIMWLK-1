//
//
//
//
// Tue Jul 18 14:24:42 EDT 2017
//

#include "JIMWLK.h"
#include "string_conv.h"

using namespace std;
using namespace blitz;

//choice of parameters g^2 \mu = 1

const int size_x=1024*2; // # steps in transverse direction
const int N_Y=100; // # steps in long direction

//const double step_x=0.05;
const double L_x=32; // transverse extent
const double step_x=L_x/size_x; // transvers step size

const double step_x2 = step_x*step_x;
const int size_x2 = size_x*size_x;

const double dY_times_alpha=0.001;

const double prefactor = sqrt(dY_times_alpha)/M_PI;
double UV = step_x*step_x*1e-5;

const int Number_of_IP=4;
double Ncoll = 6.0;
vector<double> oneOverR(Number_of_IP); 

string eventID;

boost::normal_distribution<> nd( 0.0, 1.0 / (step_x * sqrt(N_Y) ) );
boost::variate_generator< RNGType, boost::normal_distribution<> > Gaussian(rng, nd);

boost::normal_distribution<> nd_zeta( 0.0, 1.0 / step_x  );
boost::variate_generator< RNGType, boost::normal_distribution<> > Gaussian_zeta(rng, nd_zeta);

typedef blitz::TinyVector<cd,9> AdComp;
typedef blitz::TinyMatrix<cd,3,3>  colorMat;
typedef blitz::TinyMatrix<cd,8,8>  colorAdMat;
typedef Array<colorMat,2>  colorArr;
typedef Array<colorAdMat,2>  colorAdArr;

ofstream fileout,t_data;
ofstream t_data2;
ofstream t_data3;


vector<double> center_x;
vector<double> center_y;




AdComp DecomposeAdMatrix(colorAdMat U)
{
    AdComp out;
    vector<colorAdMat> F=get_F();
    for(int i=0; i<9; i++)
    {
        cd component = 1.0/3.0*trace(Product(F.at(i),U));
        out(i)=component;
    }
	//out(8) = 1.0/3.0*trace(Product(F.at(8),U));
    return out;
}

double int_to_x(int i)
{
    return i*step_x;
}

void random_center()
{
    for(int i=0;i<10;i++)
    {
        center_x.push_back(int_to_x(rand() % size_x));
        center_y.push_back(int_to_x(rand() % size_x));
    }
}

double Qs(double Y)
{
    double A0 = 0.269933;
    double A1 = 1.60715;
    double A2 = 0.430101;
    return A0*exp(A1*Y)+A2;
}


int x_to_int(double x)
{
    int out = int(x/step_x+0.5);
    return out;
}




double x2(double x)
{
    return x*x;
}


double L_d_by_2pi = L_x/(2.0*M_PI);

Array<cd,2> f_FT_Kx(void)
{
    Array<cd,2> Kx(size_x,size_x);
    Array<cd,2> FT_Kx(size_x,size_x);

    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            double x = int_to_x(i);
            double y = int_to_x(j);

            /*
            int nx_prime=0;
            int ny_prime=0;
            double distance_prime=x*x+y*y;
            for(int nx=-1;nx<2;nx++)
            for(int ny=-1;ny<2;ny++)
            {
            	double distance = (x+L_x*nx)*(x+L_x*nx)+(y+L_x*ny)*(y+L_x*ny);
            	if(distance<=distance_prime)
            	{
            		distance_prime=distance;
            		nx_prime=nx;
            		ny_prime=ny;
            	}
            }

            double x_p = x + L_x*nx_prime;
            double y_p = y + L_x*ny_prime;

            Kx(i,j) = x_p/(x_p*x_p+y_p*y_p+UV);
            */
            double denominator = ( x2(sin(0.5*x/L_d_by_2pi)) + x2(sin(0.5*y/L_d_by_2pi)) + 1e-20 ) * x2(2.0*L_d_by_2pi );
            Kx(i,j) = L_d_by_2pi * sin (x/L_d_by_2pi ) / denominator;
        }
    }

    FFTW(Kx,FT_Kx);
    return FT_Kx;
}

Array<cd,2> f_FT_Ky(void)
{
    Array<cd,2> Ky(size_x,size_x);
    Array<cd,2> FT_Ky(size_x,size_x);

    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            double x = int_to_x(i);
            double y = int_to_x(j);


            /*int nx_prime=0;
            int ny_prime=0;
            double distance_prime=x*x+y*y;
            for(int nx=-1;nx<2;nx++)
            for(int ny=-1;ny<2;ny++)
            {
            	double distance = (x+L_x*nx)*(x+L_x*nx)+(y+L_x*ny)*(y+L_x*ny);
            	if(distance<=distance_prime)
            	{
            		distance_prime=distance;
            		nx_prime=nx;
            		ny_prime=ny;
            	}
            }

            double x_p = x + L_x*nx_prime;
            double y_p = y + L_x*ny_prime;

            Ky(i,j) = y_p/(x_p*x_p+y_p*y_p+UV);
            */
            double denominator = ( x2(sin(0.5*x/L_d_by_2pi)) + x2(sin(0.5*y/L_d_by_2pi)) + 1e-20 ) * x2(2.0*L_d_by_2pi );
            Ky(i,j) = L_d_by_2pi * sin (y/L_d_by_2pi ) / denominator;

        }
    }

    FFTW(Ky,FT_Ky);
    return FT_Ky;
}



//Here 
void generate_zeta0(blitz::Array<double,2>& zeta0 )
{

    for(int i=0; i<size_x; i++)
    {
        //cout << i << "\n" << flush;
        for(int j=0; j<size_x; j++)
        {
            zeta0(i,j) = Gaussian_zeta();
        }
    }

    // Periodic boundary conditions
    //for(int i=0; i<size_x; i++) zeta0(i,0)=zeta0(i,size_x-1);
    //for(int j=0; j<size_x; j++) zeta0(0,j)=zeta0(size_x-1,j);
}



//Here 
void generate_zeta(colorArr& zeta)
{
    blitz::Array<double,2> zeta0[8];
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();

    for(int a=0; a<8; a++)
    {
        zeta0[a].resize(size_x,size_x);
        generate_zeta0(zeta0[a]);
    }


    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {

            colorMat zeta_sum;
            zeta_sum=0;

            for(int a=0; a<8; a++)
            {
                zeta_sum += zeta0[a](i,j) * lambda.at(a);
            }

            zeta(i,j) = zeta_sum;
        }
    }
}



void AdWilsonLine(colorAdArr& U, colorArr& V)
{
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();
    for(int a=0; a<8; a++)
        for(int b=0; b<8; b++)
        {
            for(int i=0; i<size_x; i++)
                for(int j=0; j<size_x; j++)
                {
                    colorMat Vx, Vxdagger, UnderTr;
                    Vx = V(i,j);
                    Vxdagger = dagger(Vx);
                    UnderTr=Product(Product(lambda.at(a),Vx), Product(lambda.at(b),Vxdagger));
                    cd Ux = 2.0*(trace(UnderTr));
                    U(i,j)(a,b) = Ux;
                }
        }
}




colorArr Integral1 (colorArr& V, colorArr& zeta_x, colorArr& zeta_y )
{
    colorArr Int(size_x,size_x);

    Array<cd,2> helper_out(size_x,size_x);

    Array<cd,2> helper_x_in(size_x,size_x);
    Array<cd,2> helper_x_out(size_x,size_x);

    Array<cd,2> helper_y_in(size_x,size_x);
    Array<cd,2> helper_y_out(size_x,size_x);



    colorArr V_zetax_Vdagger(size_x,size_x);
    colorArr V_zetay_Vdagger(size_x,size_x);



    Array<cd,2> FT_Kx(size_x,size_x);
    Array<cd,2> FT_Ky(size_x,size_x);

    FT_Kx = f_FT_Kx();
    FT_Ky = f_FT_Ky();

    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            colorMat Vdag = dagger(V(i,j));
            V_zetax_Vdagger(i,j)  =    Product(V(i,j), Product( zeta_x(i,j), Vdag) );
            V_zetay_Vdagger(i,j)  =    Product(V(i,j), Product( zeta_y(i,j), Vdag) );
        }
    }

    Array <cd,2> FT_K_V_z_Vd(size_x,size_x);

    for(int a=0; a<3; a++)
        for(int b=0; b<3; b++)
        {
            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    helper_x_in(i,j) = 	V_zetax_Vdagger(i,j)(a,b);
                    helper_y_in(i,j) = 	V_zetay_Vdagger(i,j)(a,b);
                }
            }

            FFTW(helper_x_in,helper_x_out);
            FFTW(helper_y_in,helper_y_out);

            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    //double k2 = k_x*k_x+k_y*k_y+1e-6;

                    FT_K_V_z_Vd(i,j) = (helper_x_out(i,j) * FT_Kx(i,j)+  helper_y_out(i,j) * FT_Ky(i,j));
                }
            }


            FFTW_b(FT_K_V_z_Vd,helper_out); //should be normalized by size_x^2

            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    Int(i,j)(a,b) = (helper_out(i,j)) * step_x2/double(size_x2) ;
                }
            }


        }
    return Int;
}



colorArr Integral2 (colorArr& zeta_x, colorArr& zeta_y )
{
    colorArr Int(size_x,size_x);

    Array<cd,2> helper_out(size_x,size_x);

    Array<cd,2> helper_x_in(size_x,size_x);
    Array<cd,2> helper_x_out(size_x,size_x);

    Array<cd,2> helper_y_in(size_x,size_x);
    Array<cd,2> helper_y_out(size_x,size_x);


    Array<cd,2> FT_Kx(size_x,size_x);
    Array<cd,2> FT_Ky(size_x,size_x);

    FT_Kx = f_FT_Kx();
    FT_Ky = f_FT_Ky();


    Array <cd,2> FT_K_z(size_x,size_x);

    for(int a=0; a<3; a++)
        for(int b=0; b<3; b++)
        {
            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    helper_x_in(i,j) = 	zeta_x(i,j)(a,b);
                    helper_y_in(i,j) = 	zeta_y(i,j)(a,b);
                }
            }

            FFTW(helper_x_in,helper_x_out);
            FFTW(helper_y_in,helper_y_out);

            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    FT_K_z(i,j) = (helper_x_out(i,j) * FT_Kx(i,j) +  helper_y_out(i,j) * FT_Ky(i,j)) ;
                }
            }


            FFTW_b(FT_K_z,helper_out); //should be normalized by size_x^2

            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    Int(i,j)(a,b) = (helper_out(i,j))  * step_x2/double(size_x2)  ;
                }
            }


        }
    return Int;
}




colorMat Evolution_kernel(int i, int j, colorArr& V, colorArr& Int1, colorArr& Int2)
{
    colorMat exp1;
    colorMat u_exp1;
    colorMat exp2;
    colorMat u_exp2;


    u_exp1 = -prefactor*Int1(i,j);
    u_exp2 = prefactor*Int2(i,j); //i's are in the matrix exponent!
    exp1 =  matrix_exp_an_antiH(u_exp1);
    exp2 =  matrix_exp_an_antiH(u_exp2);

    colorMat out;
    out = Product( Product( exp1,  V(i,j) ), exp2   );
    return out;
}



void evo_step(colorArr &V_prev, colorArr &V_next)
{


    colorArr zeta_x(size_x,size_x);
    colorArr zeta_y(size_x,size_x);

    generate_zeta(zeta_x);
    generate_zeta(zeta_y);


    colorArr Int1(size_x,size_x);
    colorArr Int2(size_x,size_x);
    Int1 = Integral1(V_prev, zeta_x, zeta_y);
    Int2 = Integral2(zeta_x, zeta_y);

    for(int i=0; i<size_x; i++)
        for(int j=0; j<size_x; j++)
        {

            V_next(i,j) = Evolution_kernel(i,j,V_prev, Int1, Int2);
        }
}



void rho_generator(blitz::Array<double,2>& rho)
{
    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            rho(i,j) = Gaussian(); // /sqrt(Ncoll);
        }
    }

    // Periodic boundary conditions
    //for(int i=0; i<size_x; i++) rho(i,0)=rho(i,size_x-1);
    //for(int j=0; j<size_x; j++) rho(0,j)=rho(size_x-1,j);

}





blitz::Array<complex<double>,2> A_image (blitz::Array<complex<double>,2>&  rho_image)
{
    blitz::Array<complex<double>,2> out(size_x,size_x);
    out=0.0;
    for(int i=0; i<size_x; i++)
        for(int j=0; j<size_x; j++)
            out(i,j) = - 0.5*step_x*step_x * rho_image(i,j)
                       /
                       ( cos(2.0*M_PI*i/double(size_x)) + cos(2.0*M_PI*j/double(size_x)) - 2.0 );
    out(0,0)=cd(0.0,0.0); //removing zero mode
    return out;
}


int index_x_boundary(int in)
{
    //periodic boundaries
    int i = in;
    if(i<0) i = i + size_x;
    if(i>size_x-1) i = i - size_x;
    return i;
}

TinyMatrix<cd,3,3> fV(int i, int j, vector<Array<double,2> >&  A_a)
{
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();

    blitz::TinyMatrix<cd,3,3>  V_out;

    V_out=cd(1.0,0.0),cd(0.0,0.0),cd(0.0,0.0),
    cd(0.0,0.0),cd(1.0,0.0),cd(0.0,0.0),
    cd(0.0,0.0),cd(0.0,0.0),cd(1.0,0.0);
    for(int ny=0; ny<N_Y; ny++)
    {

        blitz::TinyMatrix<cd,3,3>  in_exp;
        blitz::TinyMatrix<cd,3,3>  V;

        in_exp=cd(0.0,0.0);

        for(int a=0; a<8; a++)
            in_exp = in_exp + lambda.at(a) * A_a.at(ny*8+a)(i,j);
        in_exp=in_exp; // times i in the matrixix exponent
        V = matrix_exp_an_antiH(in_exp);
        V_out=Product(V,V_out);
    }

    return V_out;
}

blitz::Array<complex<double>,2> fft(blitz::Array<double,2>& rho)
{

    blitz::Array<complex<double>,2> image(size_x,size_x);

    blitz::Array<complex<double>,2> in(size_x,size_x);
    in = rho;
    FFTW(in,image);

    return image;
}

void IC_MV( colorArr& V )
{
    blitz::TinyMatrix<cd,3,3>  V_unit;
    V_unit=cd(1.0,0.0),cd(0.0,0.0),cd(0.0,0.0),
    cd(0.0,0.0),cd(1.0,0.0),cd(0.0,0.0),
    cd(0.0,0.0),cd(0.0,0.0),cd(1.0,0.0);

    for(int i=0; i<size_x; i=i+1)
    {
        for(int j=0; j<size_x; j=j+1)
        {
            V(i,j)=V_unit;
        }
    }

    for(int ny=0; ny<N_Y; ny++)
    {

        cerr << "target slices" << ny << "\n";

        vector<blitz::Array<double,2> >  A_a(8);


        for(int i=0; i<8; i++)
        {
            A_a.at(i).resize(size_x, size_x);
        }

        for(int i=0; i<8; i++)
        {
            blitz::Array<double,2> rho_1(size_x, size_x);
            rho_generator(rho_1);

            blitz::Array<cd,2>  image(size_x, size_x);

            image=fft(rho_1);
            rho_1.free(); 
            blitz::Array<cd,2>  tmp(size_x, size_x);
			tmp=A_image(image);
			image.free();
            blitz::Array<cd,2>  A(size_x, size_x);
            FFTW_b(tmp, A);
            tmp.free();
			A=A /  double(size_x*size_x);
            blitz::Array<double,2>  B(size_x, size_x);
            B=real(A);
            A_a.at(i)=B;
        }

        for(int i=0; i<size_x; i=i+1)
        {
            for(int j=0; j<size_x; j=j+1)
            {
                vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();


                blitz::TinyMatrix<cd,3,3>  in_exp;
                blitz::TinyMatrix<cd,3,3>  V_at_ij;

                in_exp=cd(0.0,0.0);

                for(int a=0; a<8; a++)
                {
                    in_exp = in_exp + lambda.at(a) * A_a.at(a)(i,j);
                }

                V_at_ij = matrix_exp_an_antiH(in_exp);
                V(i,j) = Product(V(i,j),V_at_ij);
            }
        }
    }

}

#include <sys/stat.h>
#include "analyzer.h"
#include "nr3.h"
#include "interp_1d.h"
#include "interp_linear.h"
#include "interp_2d.h"
#include "fourier.h"

string dname;

complex<double> su3_group_element(colorMat V,  int a)
{
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();
    colorMat P;
    P = Product(V, lambda.at(a));
    return 2.0*trace(P);
}

void output(double Y,colorArr& V_c)
{
    string fname;
    ofstream d_data;
    fname = dname+"/S_" + toString(Y) + ".dat";
    d_data.open(fname.c_str());

    for(int center=0;center<center_x.size();center++)
    {
	vector<blitz::Array<cd,2> >  comp(Number_of_IP);
    vector<blitz::Array<cd,2> >  compFT(Number_of_IP);
    vector<blitz::Array<cd,2> >  sum(Number_of_IP);
    vector<blitz::Array<cd,2> >  sum_proton(Number_of_IP);
    vector<blitz::Array<cd,2> >  sumAll(Number_of_IP);
    vector<blitz::Array<cd,2> >  pS(Number_of_IP);
    vector<blitz::Array<cd,2> >  pS_proton(Number_of_IP);
    vector<blitz::Array<cd,2> >  pSAll(Number_of_IP);

    double Qs_Y = Qs(Y);

	for (int i=0;i<Number_of_IP; i++)
	{
		comp.at(i).resize(size_x,size_x); 
		compFT.at(i).resize(size_x,size_x); 
		sum.at(i).resize(size_x,size_x);
        sum_proton.at(i).resize(size_x,size_x);
		sumAll.at(i).resize(size_x,size_x); 
		pS.at(i).resize(size_x,size_x);
        pS_proton.at(i).resize(size_x,size_x);
        pSAll.at(i).resize(size_x,size_x);
		
		comp.at(i) = 0.0;
		compFT.at(i) = 0.0;
		sum.at(i) = 0.0;
		sumAll.at(i) = 0.0;
		pSAll.at(i) = 0.0;
		pS.at(i) = 0.0;
	}

    int max_k_int = 100;
    double max_k = 15.0;
    double step_k = max_k/max_k_int;

    for(int N=0; N<Number_of_IP; N++)
    for(int a=0; a<9; a++)
    {
        for(int i=0; i<size_x; i=i+1)
        {
            for(int j=0; j<size_x; j=j+1)
            {
                double x = center_x.at(center)-int_to_x(i);
                double y = center_y.at(center)-int_to_x(j);
                
                if (x > L_x * 0.5) x = x - L_x;
                if (x <= - L_x * 0.5) x = x + L_x;

                if (y > L_x * 0.5) y = y - L_x;
                if (y <= - L_x * 0.5) y = y + L_x;
                    
                double distance2 = x*x+y*y;
                
                comp.at(N)(i,j) = su3_group_element(V_c(i,j), a) * exp(-0.5*distance2*pow(oneOverR.at(N),2)) ;
            }
        }

        FFTW(comp.at(N), compFT.at(N));

        for(int i=0; i<size_x; i=i+1)
        {
            double kx  = 2.0*M_PI*i/L_x;
            double kx_t  = 2.0/step_x*sin(kx*step_x/2.0);
            for(int j=0; j<size_x; j=j+1)
            {

                double ky  = 2.0*M_PI*j/L_x;
                double ky_t  = 2.0/step_x*sin(ky*step_x/2.0);
                double k2 = (kx_t*kx_t+ky_t*ky_t);

                sumAll.at(N)(i,j) += compFT.at(N)(i,j)*conj(compFT.at(N)(i,j));

                if(k2>pow(10.0*Qs_Y*Qs_Y/Qs(0),2))
                if(k2<1.0/step_x2/(Ncoll))
                {
                    sum.at(N)(i,j) += k2*compFT.at(N)(i,j)*conj(compFT.at(N)(i,j));
                }

                if(k2>pow(sqrt(Ncoll)*10.0*Qs_Y*Qs_Y/Qs(0),2))
                    if(k2<1.0/step_x2)
                    {
                        sum_proton.at(N)(i,j) += k2*compFT.at(N)(i,j)*conj(compFT.at(N)(i,j));
                    }
                
            }
        }
    }

    for(int i=0; i<size_x; i=i+1)
    {
        double kx  = 2.0*M_PI*i/L_x;
        double kx_t  = 2.0/step_x*sin(kx*step_x/2);
        for(int j=0; j<size_x; j=j+1)
        {
            double ky  = 2.0*M_PI*j/L_x;
            double ky_t  = 2.0/step_x*sin(ky*step_x/2);
            double k2 = (kx_t*kx_t+ky_t*ky_t);

            if( sqrt(k2) < 1.0/step_x )
            {
                fileout << Y << " " << sqrt(k2)/Qs(Y); 
				for(int N=0;N<Number_of_IP; N++)
				{
					fileout << " " << real(sumAll.at(N)(i,j))*0.5/3.0/size_x2;
				}
                fileout << " " << Qs_Y << " " << center << "\n" << flush;
            }
        }
    }

	for(int N=0;N<Number_of_IP; N++) FFTW_b(sumAll.at(N), pSAll.at(N));
	for(int N=0;N<Number_of_IP; N++) FFTW_b(sum.at(N), pS.at(N));
    for(int N=0;N<Number_of_IP; N++) FFTW_b(sum_proton.at(N), pS_proton.at(N));
    std::cout.precision(15);

    cout << Y;  // 1 
    for(int N=0;N<Number_of_IP; N++)	cout <<  " " << real(pS.at(N)(0,0))*0.5/3.0/size_x2/size_x2;
    for(int N=0;N<Number_of_IP; N++)	cout <<  " " << real(pS_proton.at(N)(0,0))*0.5/3.0/size_x2/size_x2;
    for(int N=0;N<Number_of_IP; N++)	cout <<  " " << real(pSAll.at(N)(0,0))*0.5/3.0/size_x2/size_x2; 
    cout << " " << Qs_Y << " " << center << "\n" << flush;
    }
}

double signA(int a)
{
if(a==8) return -1.0;
return 1.0;
}

void XofK (double Y, colorArr& V_c)
{
	colorAdArr Um(size_x,size_x);
   	AdWilsonLine(Um, V_c);
				
	double Qs_Y = Qs(Y);
	
	Array<cd,2> fftW(size_x, size_x);
	Array<cd,2> fieldW(size_x, size_x);
	vector<Array<cd,2> > SumW(Number_of_IP);
	vector<Array<cd,2> > bSumW(Number_of_IP);

	vector<Array<cd,2> > SumW_A(Number_of_IP);
	vector<Array<cd,2> > bSumW_A(Number_of_IP);

	vector<Array<cd,2> > SumW_p(Number_of_IP);
	vector<Array<cd,2> > bSumW_p(Number_of_IP);


    for(int center=0;center<center_x.size();center++)
    {

	for (int i=0;i<Number_of_IP; i++)
	{
		SumW.at(i).resize(size_x,size_x); 
		bSumW.at(i).resize(size_x,size_x); 
		SumW_A.at(i).resize(size_x,size_x); 
		bSumW_A.at(i).resize(size_x,size_x); 	
		SumW_p.at(i).resize(size_x,size_x); 
		bSumW_p.at(i).resize(size_x,size_x); 	
		
		SumW.at(i) = 0.0;
		bSumW.at(i) = 0.0;
		SumW_A.at(i) = 0.0;
		bSumW_A.at(i) = 0.0;
		SumW_p.at(i) = 0.0;
		bSumW_p.at(i) = 0.0;

	}

	for(int N=0;N<Number_of_IP;N++)
	for(int a=0;a<8;a++)
	for(int b=0;b<8;b++)
	{    
		fftW = cd(0.0,0.0);
    	fieldW = cd(0.0,0.0);
		for(int i=0; i<size_x; i=i+1)
    	{
        	for(int j=0; j<size_x; j=j+1)
        	{

				double x = center_x.at(center)-int_to_x(i);
                double y = center_y.at(center)-int_to_x(j);

                if (x > L_x * 0.5) x = x - L_x;
                if (x <= - L_x * 0.5) x = x + L_x;
                
                if (y > L_x * 0.5) y = y - L_x;
                if (y <= - L_x * 0.5) y = y + L_x;
                
                double distance2 = x*x+y*y;

            	fieldW(i,j) = Um(i,j)(a,b)*exp(-0.5*distance2*pow(oneOverR.at(N),2)); 
        	}
    	}

    	FFTW(fieldW, fftW);
		for(int i=0; i<size_x; i=i+1)
    	{
        	for(int j=0; j<size_x; j=j+1)
        	{
            	SumW.at(N)(i,j) = SumW.at(N)(i,j) + fftW(i,j)*conj(fftW(i,j));
        	}
    	}
	}


	double stepsize=2*2*M_PI/L_x; 
	int Nstep=int(1.0/step_x/(stepsize))+1;  
	vector<vector<double> > values(Number_of_IP); 
	vector<vector<double> > number(Number_of_IP); 
	for(int i=0;i<Number_of_IP;i++)
	{
		values.at(i).resize(Nstep); 
		number.at(i).resize(Nstep); 
		for(int j=0;j<Nstep;j++)
		{
			values.at(i).at(j)=0.0;
			number.at(i).at(j)=0.0;
		}
	}
    
    for(int i=0; i<size_x; i=i+1)
    {
        double kx  = 2.0*M_PI*i/L_x;
        double kx_t  = 2.0/step_x*sin(kx*step_x/2);
		
		for(int j=0; j<size_x; j=j+1)
        {
            double ky  = 2.0*M_PI*j/L_x;
            double ky_t  = 2.0/step_x*sin(ky*step_x/2);
            double k2 = (kx_t*kx_t+ky_t*ky_t);

            if( (sqrt(k2) < 1.0/step_x ) )
            {
				int ik = int( (sqrt(k2)) / stepsize);
				for(int N=0;N<Number_of_IP;N++)
				{
					number.at(N).at(ik)+=1.0; 
					values.at(N).at(ik)+=real(SumW.at(N)(i,j))/size_x2; 
				}


				if( (k2>pow(10.0*Qs_Y*Qs_Y/Qs(0),2)) )
				{
					for(int N=0;N<Number_of_IP;N++)
						SumW_A.at(N)(i,j) = SumW.at(N)(i,j); 
				}

				if( (k2>pow(sqrt(Ncoll)*10.0*Qs_Y*Qs_Y/Qs(0),2)) )
				{
					for(int N=0;N<Number_of_IP;N++)
						SumW_p.at(N)(i,j) = SumW.at(N)(i,j); 
				}
            }
        }
    }

	for(int ik=0;ik<Nstep;ik++)
	{
	 			t_data << Y << " " << (double(ik)+0.5)*stepsize/Qs(Y); 
				for(int N=0;N<Number_of_IP;N++)
					t_data << " " << values.at(N).at(ik)/number.at(N).at(ik); 
				t_data << " " << Qs(Y) << " "<< center<<"\n" << flush;
	}
               
	for(int N=0; N<Number_of_IP; N++) FFTW_b(SumW.at(N), bSumW.at(N));
	for(int N=0; N<Number_of_IP; N++) FFTW_b(SumW_A.at(N), bSumW_A.at(N));
	for(int N=0; N<Number_of_IP; N++) FFTW_b(SumW_p.at(N), bSumW_p.at(N));

    cout << Y;  // 1 
    for(int N=0;N<Number_of_IP; N++)	cout <<  " " << real(bSumW_A.at(N)(0,0))*0.5/3.0/size_x2/size_x2;
    for(int N=0;N<Number_of_IP; N++)	cout <<  " " << real(bSumW_p.at(N)(0,0))*0.5/3.0/size_x2/size_x2;
    for(int N=0;N<Number_of_IP; N++)	cout <<  " " << real(bSumW.at(N)(0,0))*0.5/3.0/size_x2/size_x2; 
    cout << " " << Qs_Y << " " << center << "\n" << flush;


	cout << "# Sum"; 
	for(int N=0;N<Number_of_IP;N++)
		cout << " " << real(bSumW.at(N)(0,0)) /8.0/double(size_x2)/double(size_x2); 
	cout << " " << center_x.at(center) <<  " " << center_y.at(center) << "\n" << flush;
    }
}

void oneOverR_init(double Y)
{
	double Qs_Y = Qs(Y);
	oneOverR.at(0) = 0.0; 
	oneOverR.at(1) = Qs(0)/3.0;
    oneOverR.at(2) = sqrt(Ncoll)*Qs(0)/3.0;
	oneOverR.at(3) = Qs(0)/(3.0/sqrt(2.0));
}

int main(void)
{

    AdComp out;
    vector<colorAdMat> F=get_F();

    
    cin >> eventID;
    cin >> Ncoll;
    
    srand (time(NULL));
    random_center();

    
    string name = "/efs/sipMD_" + toString(eventID) + ".dat";
    fileout.open(name.c_str());

    string fTMDname = "/efs/SIP_" + toString(eventID) + ".dat";
    t_data.open(fTMDname.c_str());

    colorArr V_c(size_x,size_x);
    colorArr V_n(size_x,size_x);



    IC_MV(V_c);
	oneOverR_init(0);
    //output(0, V_c); 
	XofK(0,V_c); 

    int i=0;
    for(double Y=0; Y<.60001; Y+=dY_times_alpha)
    {
        evo_step(V_c, V_n);
        V_c = V_n;
        //cout << Y+dY_times_alpha << "\n" << flush;
        i++;
        if(i==100)
        {
			oneOverR_init(Y+dY_times_alpha);
            //output(Y+dY_times_alpha, V_c);
			XofK(Y+dY_times_alpha,V_c); 
            i=0;
        }
    }
}
