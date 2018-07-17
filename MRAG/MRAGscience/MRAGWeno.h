#pragma once
#include "MRAGcore/MRAGCommon.h"

namespace MRAG
{
	namespace Science
	{
#define MRAGscience_WENOEPS 1.0e-6
		
		template<typename T> 
		struct Weno5
		{
			mutable T is0, is1, is2, alpha0, alpha1, alpha2, omega0, omega1, omega2;
			
			void operator()(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, T sides[2])	const
			{
				is0 = a*(a*1.33333  - b*6.33333  + c*3.66667) + b*(b*8.33333  - c*10.3333) + c*c*3.33333;//13.0/12.0 * (a-2.0*b+c)*(a-2.0*b+c) + 1.0/4.0 * (a-4.0*b+3.0*c)*(a-4.0*b+3.0*c);
				is1 = b*(b*1.33333  - c*4.33333  + d*1.66667) + c*(c*4.33333  - d*4.33333) + d*d*1.33333;//13.0/12.0 * (b-2.0*c+d)*(b-2.0*c+d) + 1.0/4.0 * (b-d)*(b-d);
				is2 = c*(c*3.33333  - d*10.3333  + e*3.66667) + d*(d*8.33333  - e*6.33333) + e*e*1.33333;//13.0/12.0 * (c-2.0*d+e)*(c-2.0*d+e) + 1.0/4.0 * (3.0*c-4.0*d+e)*(3.0*c-4.0*d+e);
				
				alpha0 = 0.1/((is0+MRAGscience_WENOEPS)*(is0+MRAGscience_WENOEPS));
				alpha1 = 0.6/((is1+MRAGscience_WENOEPS)*(is1+MRAGscience_WENOEPS));
				alpha2 = 0.3/((is2+MRAGscience_WENOEPS)*(is2+MRAGscience_WENOEPS));
				
				omega0=alpha0/(alpha0+alpha1+alpha2);
				omega1=alpha1/(alpha0+alpha1+alpha2);
				omega2=alpha2/(alpha0+alpha1+alpha2);
				
				sides[0] = omega0*(1.0/3.0*a-7.0/6.0*b+11.0/6.0*c) + omega1*(-1.0/6.0*b+5.0/6.0*c+1.0/3.0*d) + omega2*(1.0/3.0*c+5.0/6.0*d-1.0/6.0*e);
				
				is0 = d*(d*3.3333333 - e*10.333333  + f*3.6666667) + e*(e*8.3333333  - f*6.3333333) +	f*f*1.3333333;//13.0/12.0 * (f-2.0*e+d)*(f-2.0*e+d) + 1.0/4.0 * (f-4.0*e+3.0*d)*(f-4.0*e+3.0*d);
				is1 = c*(c*1.3333333 - d*4.3333333  + e*1.6666667) + d*(d*4.3333333  - e*4.3333333) +	e*e*1.3333333;//13.0/12.0 * (e-2.0*d+c)*(e-2.0*d+c) + 1.0/4.0 * (e-c)*(e-c);
				is2 = b*(b*1.3333333 - c*6.3333333  + d*3.6666667) + c*(c*8.3333333  - d*10.333333) +	d*d*3.3333333;//13.0/12.0 * (d-2.0*c+b)*(d-2.0*c+b) + 1.0/4.0 * (3.0*d-4.0*c+b)*(3.0*d-4.0*c+b);
				
				alpha0 = 0.1/((is0+MRAGscience_WENOEPS)*(is0+MRAGscience_WENOEPS));
				alpha1 = 0.6/((is1+MRAGscience_WENOEPS)*(is1+MRAGscience_WENOEPS));
				alpha2 = 0.3/((is2+MRAGscience_WENOEPS)*(is2+MRAGscience_WENOEPS));
				
				omega0=alpha0/(alpha0+alpha1+alpha2);
				omega1=alpha1/(alpha0+alpha1+alpha2);
				omega2=alpha2/(alpha0+alpha1+alpha2);		   
				
				sides[1] =  omega0*(1.0/3.0*f-7.0/6.0*e+11.0/6.0*d) + omega1*(-1.0/6.0*e+5.0/6.0*d+1.0/3.0*c) + omega2*(1.0/3.0*d+5.0/6.0*c-1.0/6.0*b);
			}
			
			void operator()(const T input[6], T sides[2]) const
			{
				is0 = input[0]*(input[0]*(4./3.)  - input[1]*(19./3.) + input[2]*(11./3.)) + input[1]*(input[1]*(25./3.)  - input[2]*(31./3.)) + input[2]*input[2]*(10./3.);//13.0/12.0 * (input[0]-2.0*input[1]+input[2])*(input[0]-2.0*input[1]+input[2]) + 1.0/4.0 * (input[0]-4.0*input[1]+3.0*input[2])*(input[0]-4.0*input[1]+3.0*input[2]);
				is1 = input[1]*(input[1]*(4./3.)  - input[2]*(13./3.)   + input[3]*(5./3.)) + input[2]*(input[2]*(13./3.)   - input[3]*(13./3.) ) + input[3]*input[3]*(4./3.);//13.0/12.0 * (input[1]-2.0*input[2]+input[3])*(input[1]-2.0*input[2]+input[3]) + 1.0/4.0 * (input[1]-input[3])*(input[1]-input[3]);
				is2 = input[2]*(input[2]*(10./3.)  - input[3]*(31./3.)  + input[4]*(11./3.)) + input[3]*(input[3]*(25./3.)  - input[4]*(19./3.)) + input[4]*input[4]*(4./3.);//13.0/12.0 * (input[2]-2.0*input[3]+input[4])*(input[2]-2.0*input[3]+input[4]) + 1.0/4.0 * (3.0*input[2]-4.0*input[3]+input[4])*(3.0*input[2]-4.0*input[3]+input[4]);
				
				alpha0 = 0.1/((is0+MRAGscience_WENOEPS)*(is0+MRAGscience_WENOEPS));
				alpha1 = 0.6/((is1+MRAGscience_WENOEPS)*(is1+MRAGscience_WENOEPS));
				alpha2 = 0.3/((is2+MRAGscience_WENOEPS)*(is2+MRAGscience_WENOEPS));
				
				omega0=alpha0/(alpha0+alpha1+alpha2);
				omega1=alpha1/(alpha0+alpha1+alpha2);
				omega2=alpha2/(alpha0+alpha1+alpha2);
				
				sides[0] = omega0*(1.0/3.0*input[0]-7.0/6.0*input[1]+11.0/6.0*input[2]) + omega1*(-1.0/6.0*input[1]+5.0/6.0*input[2]+1.0/3.0*input[3]) + omega2*(1.0/3.0*input[2]+5.0/6.0*input[3]-1.0/6.0*input[4]);
				
				is0 = input[3]*(input[3]*(10./3.) - input[4]*(31./3.)  + input[5]*(11./3.)) + input[4]*(input[4]*(25./3.)  - input[5]*(19./3.)) +	input[5]*input[5]*(4./3.);//13.0/12.0 * (input[5]-2.0*input[4]+input[3])*(input[5]-2.0*input[4]+input[3]) + 1.0/4.0 * (input[5]-4.0*input[4]+3.0*input[3])*(input[5]-4.0*input[4]+3.0*input[3]);
				is1 = input[2]*(input[2]*(4./3.) - input[3]*(13./3.)   + input[4]*(5./3.)) + input[3]*(input[3]*(13./3.)  - input[4]*(13./3.) ) +	input[4]*input[4]*(4./3.);//13.0/12.0 * (input[4]-2.0*input[3]+input[2])*(input[4]-2.0*input[3]+input[2]) + 1.0/4.0 * (input[4]-input[2])*(input[4]-input[2]);
				is2 = input[1]*(input[1]*(4./3.) - input[2]*(19./3.)  + input[3]*(11./3.)) + input[2]*(input[2]*(25./3.)  - input[3]*(31./3.)) +	input[3]*input[3]*(10./3.);//13.0/12.0 * (input[3]-2.0*input[2]+input[1])*(input[3]-2.0*input[2]+input[1]) + 1.0/4.0 * (3.0*input[3]-4.0*input[2]+input[1])*(3.0*input[3]-4.0*input[2]+input[1]);
				
				alpha0 = 0.1/((is0+MRAGscience_WENOEPS)*(is0+MRAGscience_WENOEPS));
				alpha1 = 0.6/((is1+MRAGscience_WENOEPS)*(is1+MRAGscience_WENOEPS));
				alpha2 = 0.3/((is2+MRAGscience_WENOEPS)*(is2+MRAGscience_WENOEPS));
				
				omega0=alpha0/(alpha0+alpha1+alpha2);
				omega1=alpha1/(alpha0+alpha1+alpha2);
				omega2=alpha2/(alpha0+alpha1+alpha2);		   
				
				sides[1] =  omega0*(1.0/3.0*input[5]-7.0/6.0*input[4]+11.0/6.0*input[3]) + omega1*(-1.0/6.0*input[4]+5.0/6.0*input[3]+1.0/3.0*input[2]) + omega2*(1.0/3.0*input[3]+5.0/6.0*input[2]-1.0/6.0*input[1]);
			}
		};
		
		template<typename T> 
		struct Weno5Grad
		{
			mutable T is0, is1, is2, alpha0, alpha1, alpha2, omega0, omega1, omega2;
			
			inline const T Dplus(const T& next, const T& me) const
			{ 
				return (next - me);
			}
			
			inline const T DplusDminus(const T& next, const T& me, const T& previous) const
			{ 
				return (next +	previous - me*2.0);
			}
			
			inline const T weno(const T& a, const T& b, const T& c, const T& d) const
			{
				is0 = 13.*(a-b)*(a-b) + 3.*(a-3.*b)*(a-3.*b);
				is1 = 13.*(b-c)*(b-c) + 3.*(b+c)*(b+c);
				is2 = 13.*(c-d)*(c-d) + 3.*(3.*c-d)*(3.*c-d);
				
				alpha0 = 1./((is0+MRAGscience_WENOEPS)*(is0+MRAGscience_WENOEPS));
				alpha1 = 6./((is1+MRAGscience_WENOEPS)*(is1+MRAGscience_WENOEPS));
				alpha2 = 3./((is2+MRAGscience_WENOEPS)*(is2+MRAGscience_WENOEPS));
				
				omega0 = alpha0/(alpha0 + alpha1 + alpha2);
				omega1 = alpha2/(alpha0 + alpha1 + alpha2);
				
				return 1./3*omega0*(a-2.*b+c) + 1./6*(omega1-0.5)*(b-2.*c+d);
			}
			
			void operator()(const T& z, const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const Real h, T sides[2]) const
			{
				
				const T commonterm = 1./(12.*h) * (-Dplus(b,a)+7.*Dplus(c,b)+7.*Dplus(d,c)-Dplus(e,d));
				T WENO2ndDer[2] = { 
					weno(DplusDminus(b,a,z), DplusDminus(c,b,a), DplusDminus(d,c,b), DplusDminus(e,d,c))/h,
					weno(DplusDminus(f,e,d), DplusDminus(e,d,c), DplusDminus(d,c,b), DplusDminus(c,b,a))/h
				};
				
				sides[0] = commonterm - WENO2ndDer[0];
				sides[1] = commonterm + WENO2ndDer[1];				
			}
			
			void operator()(const T input[7], const Real h, T sides[2]) const
			{
				
				const T commonterm = 1./(12.*h) * (-1.*Dplus(input[2],input[1])+7.*Dplus(input[3],input[2])+7.*Dplus(input[4],input[3])-1.*Dplus(input[5],input[4]));
				T WENO2ndDer[2] = { 
					weno(DplusDminus(input[2],input[1],input[0]), DplusDminus(input[3],input[2],input[1]), DplusDminus(input[4],input[3],input[2]), DplusDminus(input[5],input[4],input[3]))/h,
					weno(DplusDminus(input[6],input[5],input[4]), DplusDminus(input[5],input[4],input[3]), DplusDminus(input[4],input[3],input[2]), DplusDminus(input[3],input[2],input[1]))/h
				};
				
				sides[0] = commonterm - WENO2ndDer[0];
				sides[1] = commonterm + WENO2ndDer[1];				
			}
			
		};
		
		template<typename T> 
		struct Weno3
		{
			mutable T is0, is1, is2, alpha0, alpha1, alpha2, omega0, omega1;
			
			void operator()(const T& b, const T& c, const T& d, const T& e, T sides[2])	const
			{
				is0 = (c-b)*(c-b);
				is1 = (d-c)*(d-c);
				
				alpha0 = (1./3.)/((is0+MRAGscience_WENOEPS)*(is0+MRAGscience_WENOEPS));
				alpha1 = (2./3.)/((is1+MRAGscience_WENOEPS)*(is1+MRAGscience_WENOEPS));
				
				omega0 = alpha0/(alpha0+alpha1);
				omega1 = 1.-omega0;
				
				sides[0] = 1./3.*d + 5./6.*c-1./6.*b
				+ (omega0-1./3.)*(1.5*c-0.5*b)
				+ (omega1-2./3.)*(0.5*d+0.5*c);
				
				is0 = (d-e)*(d-e);
				is1 = (d-c)*(d-c);
				
				alpha0 = (1./3.)/((is0+MRAGscience_WENOEPS)*(is0+MRAGscience_WENOEPS));
				alpha1 = (2./3.)/((is1+MRAGscience_WENOEPS)*(is1+MRAGscience_WENOEPS));
				
				omega0 = alpha0/(alpha0+alpha1);
				omega1 = 1.-omega0;
				
				sides[1] = 1./3.*c + 5./6.*d -1./6.*e
				+ (omega0-1./3.)*(1.5*d-0.5*e)
				+ (omega1-2./3.)*(0.5*c+0.5*d);	
			}
			
			void operator()(const T input[4], T sides[2]) const
			{
				is0 = (input[1]-input[0])*(input[1]-input[0]);
				is1 = (input[2]-input[1])*(input[2]-input[1]);
				
				alpha0 = (1./3.)/((is0+MRAGscience_WENOEPS)*(is0+MRAGscience_WENOEPS));
				alpha1 = (2./3.)/((is1+MRAGscience_WENOEPS)*(is1+MRAGscience_WENOEPS));
				
				omega0 = alpha0/(alpha0+alpha1);
				omega1 = 1.-omega0;
				
				sides[0] = 1./3.*input[2] + 5./6.*input[1]-1./6.*input[0]
				+ (omega0-1./3.)*(1.5*input[1]-0.5*input[0])
				+ (omega1-2./3.)*(0.5*input[2]+0.5*input[1]);
				
				is0 = (input[2]-input[3])*(input[2]-input[3]);
				is1 = (input[2]-input[1])*(input[2]-input[1]);
				
				alpha0 = (1./3.)/((is0+MRAGscience_WENOEPS)*(is0+MRAGscience_WENOEPS));
				alpha1 = (2./3.)/((is1+MRAGscience_WENOEPS)*(is1+MRAGscience_WENOEPS));
				
				omega0 = alpha0/(alpha0+alpha1);
				omega1 = 1.-omega0;
				
				sides[1] = 1./3.*input[1] + 5./6.*input[2] -1./6.*input[3]
				+ (omega0-1./3.)*(1.5*input[2]-0.5*input[3])
				+ (omega1-2./3.)*(0.5*input[1]+0.5*input[2]);
			}
		};
		
		template<typename T> 
		struct Weno3Grad
		{
			mutable T is0, is1, is2, alphap, alpham, omega0, omega1;
			
			inline void weno(const T& a, const T& b, const T& c, const T& d, const T& e) const
			{
				is0 = c-2.*b+a;
				is1 = d-2.*c+b;
				is2 = e-2.*d+c;
				
				alpham = (is0*is0+MRAGscience_WENOEPS)/(is1*is1+MRAGscience_WENOEPS);
				alphap = (is2*is2+MRAGscience_WENOEPS)/(is1*is1+MRAGscience_WENOEPS);
				
				omega0 = 1./(1.+2.*alpham*alpham);
				omega1 = 1./(1.+2.*alphap*alphap);		
			}
			
			void operator()(const T& a, const T& b, const T& c, const T& d, const T& e, const Real h, T sides[2]) const
			{
				const T commonterm = d-b;
				weno(a, b, c, d, e);		
				sides[0] = 0.5*(commonterm - omega0*(d-3.*(c-b)-a))/h;
				sides[1] = 0.5*(commonterm - omega1*(e-3.*(d-c)-b))/h;
				
			}
			
			void operator()(const T input[5], const Real h, T sides[2]) const
			{
				const T commonterm = input[3]-input[1];
				weno(input[0], input[1], input[2], input[3], input[4]);		
				sides[0] = 0.5*(commonterm - omega0*(input[3]-3.*(input[2]-input[1])-input[0]))/h;
				sides[1] = 0.5*(commonterm - omega1*(input[4]-3.*(input[3]-input[2])-input[1]))/h;
				
			}
			
		};
	}
}