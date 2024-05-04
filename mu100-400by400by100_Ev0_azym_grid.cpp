#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/odeint.hpp> // odeint function definitions
#include <math.h>
#include <cmath>
#include <omp.h>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <stdexcept>
#include <boost/throw_exception.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/integrate/max_step_checker.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_const.hpp>
#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/numeric/odeint/util/copy.hpp>
#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <iostream>
#include <omp.h>
using namespace std;
using namespace boost::numeric::odeint;
typedef std::vector< double> state_type;
ofstream fp;
#define costhslices 3000
#define xslices 1
#define yslices 1
double theta=0.15;
double omega=1.0e0;
double lam=1.0e3;
double alpha=4.0e0/3.0e0;
double c=3.0e5; // speed of light in km per sec
double mup=1.0e3;
int tsize = 8*costhslices*xslices*yslices;
double costharray[costhslices];
double cosdtharray[costhslices];
double d[8*xslices*yslices];
double dc[8*xslices*yslices];
void my_observer( const state_type &x, const double t );
double normal(double x, double sig,double mean)
{
    double ans;
    //ans=1.0e0/(sqrt(2.0e0*M_PI*sig*sig));
    ans=1.0e0*exp(-(x-mean)*(x-mean)/(2.0e0*sig*sig));
    return(ans);
}

double mupt(const double t)
{
    double ans;
    //ans = mup; // constant
    ans = 3.0e0 * mup * exp(-1.0e0*c/10.0e0 * t) ; // exp, rate is 1km^-1
    return(ans);
}


void initarrays(state_type &x)
{

    // init arrays
    for(int k=0;k<costhslices;k++)
    {
        double costhmin = -1.0e0;
        double costhmax = +1.0e0;
        costharray[k] = costhmin + ((costhmax-costhmin)*(0.5e0+(double)k))*1.0e0/(double)costhslices;
        cosdtharray[k] = (costhmax-costhmin)*1.0e0/(double)costhslices;
        
    }

    for(int i=0;i<xslices;i++)
    {
        for(int j=0;j<yslices;j++)
        {
            for(int k=0;k<costhslices;k++)
            {

                x[i*yslices*costhslices*8+j*costhslices*8+k*8+0] = 1.0e0; // neutrinos
                x[i*yslices*costhslices*8+j*costhslices*8+k*8+4] = alpha; // anti-neutrinos

                x[i*yslices*costhslices*8+j*costhslices*8+k*8+1]=0.0e0;
                x[i*yslices*costhslices*8+j*costhslices*8+k*8+2]=0.0e0;
                x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]=0.0e0;   // seed?

                x[i*yslices*costhslices*8+j*costhslices*8+k*8+5]=0.0e0;
                x[i*yslices*costhslices*8+j*costhslices*8+k*8+6]=0.0e0;
                x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]=0.0e0;
            }
        }
    }

}




void my_system( const state_type &x , state_type &dxdt , const double t )
{
    //std::cout<<"t="<<t<<endl;
    for(int i=0;i<tsize;i++)
    {
	    dxdt[i]=0.0e0;
    }
    for(int i=0;i<xslices;i++)
    {
        for(int j=0;j<yslices;j++)
        {
            for(int l=0;l<8;l++)
            {
                d[i*yslices*8+j*8+l]=0.0e0;
                dc[i*yslices*8+j*8+l]=0.0e0;
            }
        }
    }
    

    int th_id;
    int nthreads;
    #pragma omp parallel private(th_id)
    {
        nthreads=omp_get_num_threads();
        th_id = omp_get_thread_num();
        //for(int i=0;i<xslices;i++)
        for(int ip=0;(ip+th_id)<xslices;ip=ip+nthreads)
        {
            int i=ip+th_id;
            for(int j=0;j<yslices;j++)
            {
                for(int k=0;k<costhslices;k++)
                {
                    for(int l=0;l<8;l++)
                    {
                        d[i*yslices*8+j*8+l]=d[i*yslices*8+j*8+l] + mupt(t)*x[i*yslices*costhslices*8+j*costhslices*8+k*8+l]*cosdtharray[k];
                        dc[i*yslices*8+j*8+l]=dc[i*yslices*8+j*8+l] + mupt(t)*costharray[k]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+l]*cosdtharray[k];
                    }
                }
            }
        }
    } //pragma loop ends here
//int th_id;
//int nthreads;
    
#pragma omp parallel private(th_id)
{
    nthreads=omp_get_num_threads();
    th_id = omp_get_thread_num();
    //std::cout<<"nthread = "<<nthreads<<"th_id = "<<th_id<<endl;
    //for(int i=0;i<xslices;i++)
    for(int ip=0;(ip+th_id)<xslices;ip=ip+nthreads)
    {
	int i=ip+th_id;
        for(int j=0;j<yslices;j++)
        {
            
	    for(int k=0;k<costhslices;k++)
            {
		
// Vacuum Hamiltonian
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]+2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]*sin(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]+-2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]*sin(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]+2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]*cos(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]+omega*(-x[i*yslices*costhslices*8+j*costhslices*8+k*8+0]*sin(2.0e0*theta) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+1]*sin(2.0e0*theta) - 2.0e0*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2]*cos(2.0e0*theta));
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]+-2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]*sin(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]+2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]*sin(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]+-2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]*cos(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]+omega*(x[i*yslices*costhslices*8+j*costhslices*8+k*8+4]*sin(2.0e0*theta) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+5]*sin(2.0e0*theta) + 2.0e0*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6]*cos(2.0e0*theta));
// Vacuum Hamiltonian
// Matter Hamiltonian
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]+0;
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]+0;
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]+-lam*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]+lam*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]+0;
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]+0;
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]+-lam*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]+lam*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6];
// Matter Hamiltonian
// Self-Interaction Hamiltonian
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]+2.0e0*d[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] - 2.0e0*d[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] - 2.0e0*d[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] + 2.0e0*d[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] - 2.0e0*costharray[k]*dc[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] + 2.0e0*costharray[k]*dc[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] + 2.0e0*costharray[k]*dc[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] - 2.0e0*costharray[k]*dc[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2];
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]+-2.0e0*d[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] + 2.0e0*d[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] + 2.0e0*d[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] - 2.0e0*d[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] + 2.0e0*costharray[k]*dc[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] - 2.0e0*costharray[k]*dc[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] - 2.0e0*costharray[k]*dc[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] + 2.0e0*costharray[k]*dc[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2];
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+0]*(d[i*yslices*8+j*8+3] - d[i*yslices*8+j*8+7] - costharray[k]*dc[i*yslices*8+j*8+3] + costharray[k]*dc[i*yslices*8+j*8+7]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+1]*(d[i*yslices*8+j*8+3] - d[i*yslices*8+j*8+7] - costharray[k]*dc[i*yslices*8+j*8+3] + costharray[k]*dc[i*yslices*8+j*8+7]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]*(d[i*yslices*8+j*8+0] - d[i*yslices*8+j*8+4] - costharray[k]*dc[i*yslices*8+j*8+0] + costharray[k]*dc[i*yslices*8+j*8+4]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]*(d[i*yslices*8+j*8+1] - d[i*yslices*8+j*8+5] - costharray[k]*dc[i*yslices*8+j*8+1] + costharray[k]*dc[i*yslices*8+j*8+5]);
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]+-x[i*yslices*costhslices*8+j*costhslices*8+k*8+0]*(d[i*yslices*8+j*8+2] - d[i*yslices*8+j*8+6] - costharray[k]*dc[i*yslices*8+j*8+2] + costharray[k]*dc[i*yslices*8+j*8+6]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+1]*(d[i*yslices*8+j*8+2] - d[i*yslices*8+j*8+6] - costharray[k]*dc[i*yslices*8+j*8+2] + costharray[k]*dc[i*yslices*8+j*8+6]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+2]*(d[i*yslices*8+j*8+0] - d[i*yslices*8+j*8+4] - costharray[k]*dc[i*yslices*8+j*8+0] + costharray[k]*dc[i*yslices*8+j*8+4]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+2]*(d[i*yslices*8+j*8+1] - d[i*yslices*8+j*8+5] - costharray[k]*dc[i*yslices*8+j*8+1] + costharray[k]*dc[i*yslices*8+j*8+5]);
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]+2.0e0*d[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] - 2.0e0*d[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] - 2.0e0*d[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] + 2.0e0*d[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] - 2.0e0*costharray[k]*dc[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] + 2.0e0*costharray[k]*dc[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] + 2.0e0*costharray[k]*dc[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] - 2.0e0*costharray[k]*dc[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6];
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]+-2.0e0*d[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] + 2.0e0*d[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] + 2.0e0*d[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] - 2.0e0*d[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] + 2.0e0*costharray[k]*dc[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] - 2.0e0*costharray[k]*dc[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] - 2.0e0*costharray[k]*dc[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] + 2.0e0*costharray[k]*dc[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6];
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+4]*(d[i*yslices*8+j*8+3] - d[i*yslices*8+j*8+7] - costharray[k]*dc[i*yslices*8+j*8+3] + costharray[k]*dc[i*yslices*8+j*8+7]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+5]*(d[i*yslices*8+j*8+3] - d[i*yslices*8+j*8+7] - costharray[k]*dc[i*yslices*8+j*8+3] + costharray[k]*dc[i*yslices*8+j*8+7]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]*(d[i*yslices*8+j*8+0] - d[i*yslices*8+j*8+4] - costharray[k]*dc[i*yslices*8+j*8+0] + costharray[k]*dc[i*yslices*8+j*8+4]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]*(d[i*yslices*8+j*8+1] - d[i*yslices*8+j*8+5] - costharray[k]*dc[i*yslices*8+j*8+1] + costharray[k]*dc[i*yslices*8+j*8+5]);
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]+-x[i*yslices*costhslices*8+j*costhslices*8+k*8+4]*(d[i*yslices*8+j*8+2] - d[i*yslices*8+j*8+6] - costharray[k]*dc[i*yslices*8+j*8+2] + costharray[k]*dc[i*yslices*8+j*8+6]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+5]*(d[i*yslices*8+j*8+2] - d[i*yslices*8+j*8+6] - costharray[k]*dc[i*yslices*8+j*8+2] + costharray[k]*dc[i*yslices*8+j*8+6]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+6]*(d[i*yslices*8+j*8+0] - d[i*yslices*8+j*8+4] - costharray[k]*dc[i*yslices*8+j*8+0] + costharray[k]*dc[i*yslices*8+j*8+4]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+6]*(d[i*yslices*8+j*8+1] - d[i*yslices*8+j*8+5] - costharray[k]*dc[i*yslices*8+j*8+1] + costharray[k]*dc[i*yslices*8+j*8+5]);
// Self-Interaction Hamiltonian
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7];
		
            }
        }
    }
}// closing parenthesis for pragma loop
}
void my_observer( const state_type &x, const double t )
{
    std::cout<<t<<" ";
    for(int i=0;i<xslices;i++)
    {
        for( int j=0;j<yslices;j++)
        {
            for(int k=0;k<costhslices;k=k+1)
            {
                for(int l=0;l<8;l++)
                {
                    std::cout<<x[i*yslices*costhslices*8+j*costhslices*8+k*8+l]<<" ";
                }
            }
        }
    }
    std::cout<<std::endl;
    
    std::string fname = __FILE__;
//    size_t pos = fname.find("cpp");
    boost::replace_all(fname, "cpp", "sum");
//    std::cout<< fname <<endl;
    ofstream fp;
    fp.open(fname,ios::app | ios::out);
    for(int i=0;i<xslices;i++)
    {
        for(int j=0;j<yslices;j++)
        {
	    for(int l=0;l<8;l++)
            {
                d[i*yslices*8+j*8+l]=0.0e0;
            }
        }
    }
    fp<<t<<" "; 
    for(int i=0;i<xslices;i++)
    {
        for(int j=0;j<yslices;j++)
        {
            for(int k=0;k<costhslices;k++)
            {
                for(int l=0;l<8;l++)
                {
                    d[i*yslices*8+j*8+l]=d[i*yslices*8+j*8+l]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+l];
                }
            }
            for(int l=0;l<8;l++)
            {
                fp<<d[i*yslices*8+j*8+l]<<" ";
            }
        }
    }
    fp<<endl;   
    fp.close();
    
    
    
        // mnr condition
        std::string ffname = "mnr.sum";
        ofstream fpp;
        fpp.open(ffname,ios::app | ios::out);
            
        double mnr_nu[xslices*yslices];
        double mnr_nu0[xslices*yslices];
    
        for(int i=0;i<xslices;i++)
        {
            for(int j=0;j<yslices;j++)
            {
                mnr_nu[i*yslices+j]=0.0e0;
                mnr_nu0[i*yslices+j]=0.0e0;
            }
        }
        double dtemp[8*xslices*yslices];
        double dctemp[8*xslices*yslices];
        for(int i=0;i<xslices;i++)
        {
            for(int j=0;j<yslices;j++)
            {
                for(int l=0;l<8;l++)
                {
                    dtemp[i*yslices*8+j*8+l]=0.0e0;
                    dctemp[i*yslices*8+j*8+l]=0.0e0;
                }
            }
        }
        

    
        for(int i=0;i<xslices;i++)
        {
            for(int j=0;j<yslices;j++)
            {
                for(int k=0;k<costhslices;k++)
                {
                    for(int l=0;l<8;l++)
                    {
                                                
                    dtemp[i*yslices*8+j*8+l]=dtemp[i*yslices*8+j*8+l] + mupt(t)*x[i*yslices*costhslices*8+j*costhslices*8+k*8+l]*cosdtharray[k];
                    dctemp[i*yslices*8+j*8+l]=dctemp[i*yslices*8+j*8+l] + mupt(t)*costharray[k]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+l]*cosdtharray[k];
                
                                        
                    }
                }
            }
        }
        
    
        fpp<<t<<" ";
        for(int i=0;i<xslices;i++)
        {
            for(int j=0;j<yslices;j++)
            {
                int kp = costhslices-1;
                // select one angle bin
                double hee    = dtemp[i*yslices*8+j*8+0] - costharray[kp]*dctemp[i*yslices*8+j*8+0];
                double hxx    = dtemp[i*yslices*8+j*8+1] - costharray[kp]*dctemp[i*yslices*8+j*8+1];
                double heebar = dtemp[i*yslices*8+j*8+4] - costharray[kp]*dctemp[i*yslices*8+j*8+4];
                double hxxbar = dtemp[i*yslices*8+j*8+5] - costharray[kp]*dctemp[i*yslices*8+j*8+5];
                
                mnr_nu[i*yslices+j] = (hee - heebar - (hxx - hxxbar));
                fpp<< mnr_nu[i*yslices+j] <<" ";
                
                
                kp = 0;
                // select other angle bin
                hee    = dtemp[i*yslices*8+j*8+0] - costharray[kp]*dctemp[i*yslices*8+j*8+0];
                hxx    = dtemp[i*yslices*8+j*8+1] - costharray[kp]*dctemp[i*yslices*8+j*8+1];
                heebar = dtemp[i*yslices*8+j*8+4] - costharray[kp]*dctemp[i*yslices*8+j*8+4];
                hxxbar = dtemp[i*yslices*8+j*8+5] - costharray[kp]*dctemp[i*yslices*8+j*8+5];
                
                mnr_nu0[i*yslices+j] = (hee - heebar - (hxx - hxxbar));
                fpp<< mnr_nu0[i*yslices+j] <<" ";
                
                    
            }
        }
//        fpp<<endl;
    
        // here extra bit
        double mnr_lam[xslices*yslices];
        for(int i=0;i<xslices;i++)
        {
            for(int j=0;j<yslices;j++)
            {
                mnr_lam[i*yslices+j]=0.0e0;
            }
        }
        for(int i=0;i<xslices;i++)
        {
            for(int j=0;j<yslices;j++)
            {
                                
                mnr_lam[i*yslices+j] =  lam;

                fpp<< mnr_lam[i*yslices+j] <<" ";
                            
            }
        }
    
        double mnr_vac[xslices*yslices];
        for(int i=0;i<xslices;i++)
        {
            for(int j=0;j<yslices;j++)
            {
                mnr_vac[i*yslices+j]=0.0e0;
            }
        }
        for(int i=0;i<xslices;i++)
        {
            for(int j=0;j<yslices;j++)
            {
                        
                mnr_vac[i*yslices+j] =  omega; // vac osc freq

                fpp<< mnr_vac[i*yslices+j] <<" ";
                            
            }
        }
    
        fpp<<endl;
        // ends extra bit
    
        fpp.close();
        // end mnr condition
    
    
    
        // moments extraction (spherical harmonics)
        std::string fffname = "moments.sum";
        ofstream fppp;
        fppp.open(fffname,ios::app | ios::out);

        double N0[8*xslices*yslices];
        double N1[8*xslices*yslices];
        double N2[8*xslices*yslices];
        double N3[8*xslices*yslices];
        double N4[8*xslices*yslices];
        double N5[8*xslices*yslices];
    
        for(int i=0;i<xslices;i++)
        {
            for(int j=0;j<yslices;j++)
            {
                for(int l=0;l<8;l++)
                {
                    N0[i*yslices*8+j*8+l]=0.0e0;
                    N1[i*yslices*8+j*8+l]=0.0e0;
                    N2[i*yslices*8+j*8+l]=0.0e0;
                    N3[i*yslices*8+j*8+l]=0.0e0;
                    N4[i*yslices*8+j*8+l]=0.0e0;
                    N5[i*yslices*8+j*8+l]=0.0e0;
                }
            }
        }
        
        fppp<<t<<" ";
        for(int i=0;i<xslices;i++)
        {
            for(int j=0;j<yslices;j++)
            {
                for(int k=0;k<costhslices;k++)
                {
                    for(int l=0;l<8;l++)
                    {
                    
                        double v1 = costharray[k];
                        double v2 = costharray[k] * costharray[k];
                        double v3 = costharray[k] * costharray[k] * costharray[k];
                        double v4 = costharray[k] * costharray[k] * costharray[k] * costharray[k];
                        double v5 = costharray[k] * costharray[k] * costharray[k] * costharray[k] * costharray[k];
                        
                        /// moments
                        //N0[i*yslices*8+j*8+l] = N0[i*yslices*8+j*8+l] + 1.0 * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        //N1[i*yslices*8+j*8+l] = N1[i*yslices*8+j*8+l] + v1 * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        //N2[i*yslices*8+j*8+l] = N2[i*yslices*8+j*8+l] + v2 * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        //N3[i*yslices*8+j*8+l] = N3[i*yslices*8+j*8+l] + v3 * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        //N4[i*yslices*8+j*8+l] = N4[i*yslices*8+j*8+l] + v4 * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        //N5[i*yslices*8+j*8+l] = N5[i*yslices*8+j*8+l] + v5 * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        
                        // legendre polynomials in v
                        N0[i*yslices*8+j*8+l] = N0[i*yslices*8+j*8+l] + 1.0 * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        N1[i*yslices*8+j*8+l] = N1[i*yslices*8+j*8+l] + v1 * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        N2[i*yslices*8+j*8+l] = N2[i*yslices*8+j*8+l] + (1.0/2.0)*(3.00*v2 - 1.00) * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        N3[i*yslices*8+j*8+l] = N3[i*yslices*8+j*8+l] + (1.0/2.0)*(5.00*v3 - 3.00*v1) * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        N4[i*yslices*8+j*8+l] = N4[i*yslices*8+j*8+l] + (1.0/8.0)*(35.0*v4 - 30.0*v2 + 3.0) * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        N5[i*yslices*8+j*8+l] = N5[i*yslices*8+j*8+l] + (1.0/8.0)*(63.0*v5 - 70.0*v3 + 15.0*v1) * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        
                        // spherical harmonics (only m=0 modes)
                        //N0[i*yslices*8+j*8+l] = N0[i*yslices*8+j*8+l] + sqrt(1.0/(4.0*M_PI)) * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        //N1[i*yslices*8+j*8+l] = N1[i*yslices*8+j*8+l] + sqrt(3.0/(4.0*M_PI)) * v1 * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        //N2[i*yslices*8+j*8+l] = N2[i*yslices*8+j*8+l] + (1.0/2.0) * sqrt(5.0/(4.0*M_PI)) * (3.00*v2 - 1.00) * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        //N3[i*yslices*8+j*8+l] = N3[i*yslices*8+j*8+l] + (1.0/2.0) * sqrt(7.0/(4.0*M_PI)) * (5.00*v3 - 3.00*v1) * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        //N4[i*yslices*8+j*8+l] = N4[i*yslices*8+j*8+l] + (3.0/16.0) * sqrt(1.0/M_PI) * (35.0*v4 - 30.0*v2 + 3.0) * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                        //N5[i*yslices*8+j*8+l] = N5[i*yslices*8+j*8+l] + (1.0/16.0) * sqrt(11.0/M_PI) * (63.0*v5 - 70.0*v3 + 15.0*v1) * x[i*yslices*costhslices*8+j*costhslices*8+k*8+l] * cosdtharray[k];
                                        
                    }
                }
                
                // off-diagonal ex (norm computed for me)
                fppp << sqrt( N0[i*yslices*8+j*8+2]*N0[i*yslices*8+j*8+2] + N0[i*yslices*8+j*8+3]*N0[i*yslices*8+j*8+3] ) << " "; // rhoexN0
                fppp << sqrt( N1[i*yslices*8+j*8+2]*N1[i*yslices*8+j*8+2] + N1[i*yslices*8+j*8+3]*N1[i*yslices*8+j*8+3] ) << " "; // rhoexN1
                fppp << sqrt( N2[i*yslices*8+j*8+2]*N2[i*yslices*8+j*8+2] + N2[i*yslices*8+j*8+3]*N2[i*yslices*8+j*8+3] ) << " "; // rhoexN2
                fppp << sqrt( N3[i*yslices*8+j*8+2]*N3[i*yslices*8+j*8+2] + N3[i*yslices*8+j*8+3]*N3[i*yslices*8+j*8+3] ) << " "; // rhoexN3
                fppp << sqrt( N4[i*yslices*8+j*8+2]*N4[i*yslices*8+j*8+2] + N4[i*yslices*8+j*8+3]*N4[i*yslices*8+j*8+3] ) << " "; // rhoexN4
                fppp << sqrt( N5[i*yslices*8+j*8+2]*N5[i*yslices*8+j*8+2] + N5[i*yslices*8+j*8+3]*N5[i*yslices*8+j*8+3] ) << " "; // rhoexN5
                
                // diagonal xx
                fppp << N0[i*yslices*8+j*8+1]  << " "; // rhoxxN0
                fppp << N1[i*yslices*8+j*8+1]  << " "; // rhoxxN1
                fppp << N2[i*yslices*8+j*8+1]  << " "; // rhoxxN2
                fppp << N3[i*yslices*8+j*8+1]  << " "; // rhoxxN3
                fppp << N4[i*yslices*8+j*8+1]  << " "; // rhoxxN4
                fppp << N5[i*yslices*8+j*8+1]  << " "; // rhoxxN5
                
                // extra bit
                // REAL off-diagonal ex
                fppp << N0[i*yslices*8+j*8+2]  << " "; // RErhoexN0
                fppp << N1[i*yslices*8+j*8+2]  << " "; // RErhoexN1
                fppp << N2[i*yslices*8+j*8+2]  << " "; // RErhoexN2
                fppp << N3[i*yslices*8+j*8+2]  << " "; // RErhoexN3
                fppp << N4[i*yslices*8+j*8+2]  << " "; // RErhoexN4
                fppp << N5[i*yslices*8+j*8+2]  << " "; // RErhoexN5
                
                // IMAG off-diagonal ex
                fppp << N0[i*yslices*8+j*8+3]  << " "; // IMrhoexN0
                fppp << N1[i*yslices*8+j*8+3]  << " "; // IMrhoexN1
                fppp << N2[i*yslices*8+j*8+3]  << " "; // IMrhoexN2
                fppp << N3[i*yslices*8+j*8+3]  << " "; // IMrhoexN3
                fppp << N4[i*yslices*8+j*8+3]  << " "; // IMrhoexN4
                fppp << N5[i*yslices*8+j*8+3]  << " "; // IMrhoexN5
                
                // diagonal ee
                fppp << N0[i*yslices*8+j*8+0]  << " "; // rhoeeN0
                fppp << N1[i*yslices*8+j*8+0]  << " "; // rhoeeN1
                fppp << N2[i*yslices*8+j*8+0]  << " "; // rhoeeN2
                fppp << N3[i*yslices*8+j*8+0]  << " "; // rhoeeN3
                fppp << N4[i*yslices*8+j*8+0]  << " "; // rhoeeN4
                fppp << N5[i*yslices*8+j*8+0]  << " "; // rhoeeN5
                
                // SAME BUT FOR BARS
                
                // diagonal xxbar
                fppp << N0[i*yslices*8+j*8+5]  << " "; // rhoxxbarN0
                fppp << N1[i*yslices*8+j*8+5]  << " "; // rhoxxbarN1
                fppp << N2[i*yslices*8+j*8+5]  << " "; // rhoxxbarN2
                fppp << N3[i*yslices*8+j*8+5]  << " "; // rhoxxbarN3
                fppp << N4[i*yslices*8+j*8+5]  << " "; // rhoxxbarN4
                fppp << N5[i*yslices*8+j*8+5]  << " "; // rhoxxbarN5
                
                // extra bit
                // REAL off-diagonal ex
                fppp << N0[i*yslices*8+j*8+6]  << " "; // RErhoexbarN0
                fppp << N1[i*yslices*8+j*8+6]  << " "; // RErhoexbarN1
                fppp << N2[i*yslices*8+j*8+6]  << " "; // RErhoexbarN2
                fppp << N3[i*yslices*8+j*8+6]  << " "; // RErhoexbarN3
                fppp << N4[i*yslices*8+j*8+6]  << " "; // RErhoexbarN4
                fppp << N5[i*yslices*8+j*8+6]  << " "; // RErhoexbarN5
                
                // IMAG off-diagonal ex
                fppp << N0[i*yslices*8+j*8+7]  << " "; // IMrhoexbarN0
                fppp << N1[i*yslices*8+j*8+7]  << " "; // IMrhoexbarN1
                fppp << N2[i*yslices*8+j*8+7]  << " "; // IMrhoexbarN2
                fppp << N3[i*yslices*8+j*8+7]  << " "; // IMrhoexbarN3
                fppp << N4[i*yslices*8+j*8+7]  << " "; // IMrhoexbarN4
                fppp << N5[i*yslices*8+j*8+7]  << " "; // IMrhoexbarN5
                
                // diagonal ee
                fppp << N0[i*yslices*8+j*8+4]  << " "; // rhoeebarN0
                fppp << N1[i*yslices*8+j*8+4]  << " "; // rhoeebarN1
                fppp << N2[i*yslices*8+j*8+4]  << " "; // rhoeebarN2
                fppp << N3[i*yslices*8+j*8+4]  << " "; // rhoeebarN3
                fppp << N4[i*yslices*8+j*8+4]  << " "; // rhoeebarN4
                fppp << N5[i*yslices*8+j*8+4]  << " "; // rhoeebarN5
                    
                
            }
        }
        fppp<<endl;

        fppp.close();
        // end moments
        
}


int main(int argc, char* argv[])
{
    //-----------
    //std::cout<< __FILE__ <<endl;
    std::string fname = __FILE__;
//    size_t pos = fname.find("cpp");
    boost::replace_all(fname, "cpp", "sum");
    //std::cout<< fname <<endl;
    //ofstream fp;
    fp.open(fname);
    //-------------
    std::cout.setf ( std::ios::scientific, std::ios::floatfield );
    std::cout.precision(15);
    int size=tsize;
    state_type x0(size); // Initial condition, vector of 2 elements (position and velocity)
    double err_abs = 1.0e-16;
    double err_rel = 1.0e-12;
    initarrays(x0);
    // Integration parameters
    double t0 = 0.0e0; // time in seconds
    int tsteps = 1000; 
    double t1 = 1.0e-04; //1.0e-07; // time in seconds
    double dt = (t1-t0)/((double)tsteps);
    typedef runge_kutta_fehlberg78< state_type > solver;
    my_observer(x0,t0);
    for(int ts=0;ts<tsteps;ts++)
    {
        integrate_adaptive( make_controlled( err_abs , err_rel , solver() ), my_system, x0 , t0+ts*dt , t0+(1+ts)*dt , dt*1.0e-1,null_observer());
        my_observer(x0,t0+(1+ts)*dt);
    }
    fp.close();

}
