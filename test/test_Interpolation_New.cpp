#include <gtest/gtest.h>

#include "Interpolation.h"
#include "ODESolver.h"



// Constants : range of NMAX, for every ACCURACY tests requiring order check in the cpp file

const int LOWNMAX = 5;
const int HIGHNMAX = 10;




// Functions claim

void Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(int given_PMAX, int given_Lagr_PMAX, int smallNMAX, int largeNMAX, int DIM, int VEC_NUM,
    std::vector<std::function<double(double, int)>> init_func_current, 
    std::vector<std::function<double(double, int)>> init_func_next,
    std::function<double(std::vector<double> /*u_now*/, std::vector<double> /*u_next*/, int, int)> func_flux, bool thesparse);

void Interpolation_New_Test__Lagr_General_Check_Order_System_Separable__New(int given_PMAX, int given_Lagr_PMAX, int smallNMAX, int largeNMAX, int DIM, int VEC_NUM,
    std::vector<std::function<double(double, int)>> init_func_current, 
    std::vector<std::function<double(double, int)>> init_func_next,
    std::function<double(std::vector<double> /* u_now */, std::vector<double> /* u_next */, int, int)> func_flux, bool thesparse);

void Substract_fu_norm(DGSolution & dg_solu1, DGSolution & dg_solu2, double & val);

double mypow(double x,int p) {
    if (p == 0) return 1.; 
    else {
        double ret = 1.;
        for (auto i = 1; i <= p; i++) ret=ret*x;
        return ret;
    }
}







// Test of 1D system

TEST(Interpolation_New_Test, Lagr_General_Check_Order_1D_System)
{
    auto init_func_1_0 = [](double x, int d)->double { assert(d==0); return std::pow((x+1.)*0.5,10); };
    auto init_func_1_1 = [](double x, int d)->double { assert(d==0); return std::sin(x+1.)/(x+2.); };
    auto init_func_1_2 = [](double x, int d)->double { assert(d==0); return std::log2(std::cos(Const::PI*x) + 2.); };
    std::vector<std::function<double(double, int)>> init_func_next = {init_func_1_0, init_func_1_1, init_func_1_2};

    auto init_func_2_0 = [](double x, int d)->double { assert(d==0); return std::sin(x*2.*Const::PI); };
    auto init_func_2_1 = [](double x, int d)->double { assert(d==0); return std::log10(x*x*x+1.)*(x-2.); };
    auto init_func_2_2 = [](double x, int d)->double { assert(d==0); return ((((-2.*x + 1.)*x + 3.)*x - 4.)*x+5.)*x/100.; };
    std::vector<std::function<double(double, int)>> init_func_current = {init_func_2_0, init_func_2_1, init_func_2_2};

    auto func_flux = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d==0 && i>=0 && i<=2);
        switch (i)
        {
        case 0:
            return std::exp( u_now[0]*u_now[0] )*std::cos( u_next[2]*u_now[1]*u_next[2] );
            break;
        case 1:
            return u_next[2]/(1.5 + std::sin(u_now[2] - 2.*u_next[1] + u_next[0]));
            break;
        case 2:
            return (u_next[2]-1.) * (u_next[1] + u_now[2]) /( u_now[0] * u_now[0] + 1. );
            break;
        default:
            break;
        }
    };

    for (auto gPMAX = 1; gPMAX <= 3; gPMAX++) 
        Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(gPMAX,gPMAX,LOWNMAX,HIGHNMAX,1,3,init_func_current,init_func_next,func_flux,true);
}
















// Test of 2D Scalar case

void Lagr_General_Check_Order_2D_Scalar_PolyInit_Template(int gPMAX)
{
    auto init_func_1 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return x*x/3.+1.-2.*x;
        else return 2*x*x-1;
    };
    std::vector<std::function<double(double, int)>> init_func_next{init_func_1};
    auto init_func_2 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return x*x+x;
        else return 2.-x*x+x/3.;
    };
    std::vector<std::function<double(double, int)>> init_func_current{init_func_2};
    auto func_flux = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d>=0 && d<=1 && i==0);
        if (d == 0)
            return std::exp(u_next[0] * u_now[0]) / ( std::exp(u_now[0]) + 2.*std::exp(u_next[0]) ) ;
        else 
            return std::log10( 11. + std::sin( u_now[0]*u_next[0] ) );
    };

    Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(gPMAX,gPMAX+1,LOWNMAX,HIGHNMAX,2,1,init_func_current,init_func_next,func_flux,true);
}

TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_PolyInit_P2) { Lagr_General_Check_Order_2D_Scalar_PolyInit_Template(2); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_PolyInit_P3) { Lagr_General_Check_Order_2D_Scalar_PolyInit_Template(3); }



void Lagr_General_Check_Order_2D_Scalar_PolyFlux_Template(int gPMAX)
{
    auto init_func_1 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::log2(2.*x+1.);
        else return std::exp(2*x)/5.;
    };
    std::vector<std::function<double(double, int)>> init_func_next{init_func_1};
    auto init_func_2 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return (1+2.*x)/(1.+x*x);
        else return (1.+x*x)/(1+2.*x);
    };
    std::vector<std::function<double(double, int)>> init_func_current{init_func_2};
    auto func_flux = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d>=0 && d<=1 && i==0);
        if (d == 0)
            return std::pow(u_next[0] - u_now[0],2) ;
        else 
            return (2.+u_now[0])*(1-2.*u_next[0]);
    };

    Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(gPMAX,gPMAX+1,LOWNMAX,HIGHNMAX,2,1,init_func_current,init_func_next,func_flux,true);
}

TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_PolyFlux_P1) {  Lagr_General_Check_Order_2D_Scalar_PolyFlux_Template(1);  }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_PolyFlux_P2) {  Lagr_General_Check_Order_2D_Scalar_PolyFlux_Template(2);  }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_PolyFlux_P3) {  Lagr_General_Check_Order_2D_Scalar_PolyFlux_Template(3);  }




void Lagr_General_Check_Order_2D_Scalar_PolyFlux_2_Template(int gPMAX)
{
    auto nonlin = [](double x)->double { return -std::log10(1.-x*0.8); };
    auto init_func_1 = [&](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return mypow(2.*x+1.,5)/100. + nonlin(x);
        else return (2.+x*x*(x*x-x+1.))*x*(1-x)/5. - nonlin(x);
    };
    std::vector<std::function<double(double, int)>> init_func_next{init_func_1};
    auto init_func_2 = [&](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return (x*(x*(1+2.*x)-4)+6)*(2.+x*x)/20. - nonlin(x);
        else return ((((1-2.*x)*x+3.)*x-4)*x+5)/20. + mypow(x+1,6)/35. + nonlin(x);
    };
    std::vector<std::function<double(double, int)>> init_func_current{init_func_2};
    auto func_flux = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d>=0 && d<=1 && i==0);
        if (d == 0)
            return mypow(u_next[0] - u_now[0],5) + mypow(u_next[0] , 4) ;
        else 
            return u_now[0]*(mypow(u_now[0] + u_next[0],4)-2.*u_next[0]);
    };

    Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(gPMAX,gPMAX+1,LOWNMAX,HIGHNMAX,2,1,init_func_current,init_func_next,func_flux,true);
}

TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_PolyFlux_2_P1) { Lagr_General_Check_Order_2D_Scalar_PolyFlux_2_Template(1); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_PolyFlux_2_P2) { Lagr_General_Check_Order_2D_Scalar_PolyFlux_2_Template(2); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_PolyFlux_2_P3) { Lagr_General_Check_Order_2D_Scalar_PolyFlux_2_Template(3); }




void Lagr_General_Check_Order_2D_Scalar_1_Template(int gPMAX)
{
    auto init_func_1 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::exp(std::cos(x+2));
        else return std::sin(x*x-2*x);
    };
    std::vector<std::function<double(double, int)>> init_func_next{init_func_1};
    auto init_func_2 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::exp(x*x+x);
        else return std::cos(2.-x*x);
    };
    std::vector<std::function<double(double, int)>> init_func_current{init_func_2};
    auto func_flux = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d>=0 && d<=1 && i==0);
        if (d == 0)
            return std::sin(u_next[0] - u_now[0] * u_next[0]) ;
        else 
            return std::cos(u_now[0]*u_next[0]);
    };

    Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(gPMAX,gPMAX+1,LOWNMAX,HIGHNMAX,2,1,init_func_current,init_func_next,func_flux,true);
}

TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_1_P1) { Lagr_General_Check_Order_2D_Scalar_1_Template(1); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_1_P2) { Lagr_General_Check_Order_2D_Scalar_1_Template(2); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_1_P3) { Lagr_General_Check_Order_2D_Scalar_1_Template(3); }




void Lagr_General_Check_Order_2D_Scalar_2_Template(int gPMAX)
{
    auto init_func_1 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::exp(std::cos(x+2))/(1.+x);
        else return std::sin(std::pow(x+1. , 0.5));
    };
    std::vector<std::function<double(double, int)>> init_func_next{init_func_1};
    auto init_func_2 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::exp( x*x*std::pow(x+1.,0.5) );
        else return (2.+x)/(1.+x*x*x);
    };
    std::vector<std::function<double(double, int)>> init_func_current{init_func_2};
    auto func_flux = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d>=0 && d<=1 && i==0);
        if (d == 0)
            return std::sin(u_next[0] * u_now[0])/(2. + std::pow(u_next[0] - 2.*u_now[0], 2.));
        else 
            return std::log10(2. + std::pow(2.*u_next[0] + u_now[0], 2.));
    };
    
    Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(gPMAX,gPMAX+1,LOWNMAX,HIGHNMAX,2,1,init_func_current,init_func_next,func_flux,true);
}

TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_2_P1) { Lagr_General_Check_Order_2D_Scalar_2_Template(1); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_2_P2) { Lagr_General_Check_Order_2D_Scalar_2_Template(2); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_Scalar_2_P3) { Lagr_General_Check_Order_2D_Scalar_2_Template(3); }





















// Test of 2D System

void Lagr_General_Check_Order_2D_System_Template(int gPMAX)
{
    auto init_func_11 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::exp(std::cos(x+2)/(1.+x));
        else return std::sin(x*x+2*x);
    };
    auto init_func_12 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return (1.+x);
        else return 1.-2./(x*x+2*x+2.);
    };
    std::vector<std::function<double(double, int)>> init_func_next{init_func_11,init_func_12};
    auto init_func_21 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::exp(-x*x/2.);
        else return (1.-x)/(2.+x*x);
    };
    auto init_func_22 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::sin(x);
        else return x*std::cos(x*2)-1.;
    };
    std::vector<std::function<double(double, int)>> init_func_current{init_func_21,init_func_22};
    auto func_flux = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d>=0 && d<=1 && i>=0 && i<=1);
        if (i == 0) 
        {
            if (d == 0)
                return std::sin(u_next[0] * u_now[1] * u_now[1])/(1. + std::pow(u_next[1] - 2.*u_now[0], 2.));
            else 
                return std::log2(1. + 0.25*std::pow(u_next[0] - 2. * u_next[1] + u_now[0], 2.));
        }
        else
        {
            if (d == 0)
                return std::exp( - u_next[1] * u_now[0] * u_now[0]);
            else 
                return std::sin(1. + std::pow(u_next[1] - u_next[0] + 2.*u_now[1], 2.));           
        }
    };

    Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(gPMAX,gPMAX+1,LOWNMAX,HIGHNMAX,2,2,init_func_current,init_func_next,func_flux,true);
}

TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_System_P1) { Lagr_General_Check_Order_2D_System_Template(1); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_System_P2) { Lagr_General_Check_Order_2D_System_Template(2); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_System_P3) { Lagr_General_Check_Order_2D_System_Template(3); }




void Lagr_General_Check_Order_2D_System_PolyFlux_Template(int gPMAX)
{
    auto init_func_11 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::cos(x*8);
        else return std::sin(x*x*10);
    };
    auto init_func_12 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::log2(1.+mypow(x,4));
        else return std::exp(-x*0.5);
    };
    std::vector<std::function<double(double, int)>> init_func_next{init_func_11,init_func_12};
    auto init_func_21 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::exp(mypow(x+1,2)/4.);
        else return std::log2(1.+x*x);
    };
    auto init_func_22 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::sin(x*10);
        else return std::cos(x*x*8);
    };
    std::vector<std::function<double(double, int)>> init_func_current{init_func_21,init_func_22};
    auto func_flux = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d>=0 && d<=1 && i>=0 && i<=1);
        if (i == 0) 
        {
            if (d == 0)
                return mypow(u_now[0]+u_next[1],3);
            else 
                return (u_now[1]-2.*u_next[0]+u_next[1])*u_now[0]*u_now[0];
        }
        else
        {
            if (d == 0)
                return mypow(u_next[0]+u_now[1],3);
            else 
                return (u_next[0]-2.*u_now[1]+u_now[0])*u_next[1]*u_next[1];           
        }
    };

    Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(gPMAX,gPMAX+1,LOWNMAX,HIGHNMAX,2,2,init_func_current,init_func_next,func_flux,true);
}

TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_System_PolyFlux_P1) { Lagr_General_Check_Order_2D_System_PolyFlux_Template(1); } 
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_System_PolyFlux_P2) { Lagr_General_Check_Order_2D_System_PolyFlux_Template(2); } 
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_System_PolyFlux_P3) { Lagr_General_Check_Order_2D_System_PolyFlux_Template(3); } 




void Lagr_General_Check_Order_2D_System_3Unknowns_Template(int gPMAX)
{
    auto init_func_11 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return -2./(1.+x);
        else return std::sin(x*x+2*x);
    };
    auto init_func_12 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::exp(std::cos(x+2));
        else return std::cos(x*x/3.)-std::cos(x*x/2.);
    };
    auto init_func_13 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return 2./(x*x+2*x+2.);
        else return std::cos(std::exp(-x));
    };
    std::vector<std::function<double(double, int)>> init_func_next{init_func_11,init_func_12,init_func_13};
    auto init_func_21 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::exp(x*x);
        else return (1.-x*x)/(3.+x);
    };
    auto init_func_22 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::sin(std::sqrt(x+1.));
        else return std::log2(x+1.);
    };
    auto init_func_23 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::log10(mypow(x,4)+1);
        else return std::sin(1./(1.+x));
    };
    std::vector<std::function<double(double, int)>> init_func_current{init_func_21,init_func_22,init_func_23};
    auto func_flux = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d>=0 && d<=1 && i>=0 && i<=2);
        if (i == 0) 
        {
            if (d == 0)
                return std::sin(u_next[0] * u_now[1] * u_now[1]);
            else
                return std::exp( - std::pow(u_next[1] - 2.*u_now[0], 2.) * 0.2 );
        }
        else if (i == 1)
        {
            if (d == 0)
                return std::log2(1. + std::exp(u_next[1] * u_now[0]));
            else 
                return std::cos(1. + 0.06*std::pow(u_next[1] - u_next[0] + 2.*u_now[1], 2.));           
        }
        else
        {
            if (d == 0)
                return std::exp(std::cos( u_next[1] - 2.*u_now[1] + u_now[0] ));
            else
                return std::log10(1. + 0.1*mypow(u_next[0] - u_now[0], 2.));
        }
    };

    Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(gPMAX,gPMAX+1,LOWNMAX,HIGHNMAX,2,3,init_func_current,init_func_next,func_flux,true);
}

TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_System_3Unknowns_P1) { Lagr_General_Check_Order_2D_System_3Unknowns_Template(1); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_System_3Unknowns_P2) { Lagr_General_Check_Order_2D_System_3Unknowns_Template(2); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_2D_System_3Unknowns_P3) { Lagr_General_Check_Order_2D_System_3Unknowns_Template(3); }

























// Test of 3D system

void Lagr_General_Check_Order_3D_System_Template(int gPMAX)
{
    auto init_func_11 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=2);
        if (d == 0) return -2./(1.+x);
        else if (d == 1) return std::cos(std::exp(x));
        else return std::sin(x*x+2*x);
    };
    auto init_func_12 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=2);
        if (d == 0) return std::exp(std::cos(x+2));
        else if (d == 1) return 2./(x*x+2*x+2.);
        else return std::cos(x*x/3.)-std::cos(x*x/2.);
    };
    std::vector<std::function<double(double, int)>> init_func_next{init_func_11,init_func_12};
    auto init_func_21 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=2);
        if (d == 0) return std::exp(x*x/5.);
        else if (d == 1) return std::sin(1./(1.+x));
        else return (1.-x)/(3.+x*x);
    };
    auto init_func_22 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=2);
        if (d == 0) return std::sin(std::sqrt(x+1.));
        else if (d == 1) return std::log10(mypow(x*0.5,4)+1);
        else return std::log2(x+1.)-0.5;
    };
    std::vector<std::function<double(double, int)>> init_func_current{init_func_21,init_func_22};
    auto func_flux = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d>=0 && d<=2 && i>=0 && i<=1);
        if (i == 0) 
        {
            if (d == 0)
                return std::sin(u_next[0] * u_now[1] * u_now[1]);
            else if (d == 1)
                return std::log10(1. + 0.1*mypow(u_next[0] - u_now[0], 2.));
            else
                return std::exp( - std::pow(u_next[1] - 2.*u_now[0], 2.) * 0.2 );
        }
        else
        {
            if (d == 0)
                return std::exp( - u_next[1] * u_now[0] * u_now[0]);
            else if (d == 1)
                return std::cos( u_next[1] - 2.*u_now[1] );
            else 
                return std::sin(1. + 0.06*std::pow(u_next[1] - u_next[0] + 2.*u_now[1], 2.));           
        }
    };

    Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(gPMAX,gPMAX+1,LOWNMAX,HIGHNMAX,3,2,init_func_current,init_func_next,func_flux,true);
}

TEST(Interpolation_New_Test, Lagr_General_Check_Order_3D_System_P1) { Lagr_General_Check_Order_3D_System_Template(1); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_3D_System_P2) { Lagr_General_Check_Order_3D_System_Template(2); }
TEST(Interpolation_New_Test, Lagr_General_Check_Order_3D_System_P3) { Lagr_General_Check_Order_3D_System_Template(3); }






















/// Lagr_Consistency: compare the 2 versions of LagrInterpolation::nonlinear_Lagr; 0 difference is expected
TEST(Interpolation_New_Test, Lagr_Consistency)
{

    // define static variables
    int DIM = 2;
    int VEC_NUM = 2;
    int NMAX = 5;
    AlptBasis::PMAX = 2;

	LagrBasis::PMAX = AlptBasis::PMAX;
	LagrBasis::msh_case = 1;

	HermBasis::PMAX = 3;    // Hermite is not tested here though
	HermBasis::msh_case = 1;

    int common_num_gauss_pt = 1+(int)(0.5*(AlptBasis::PMAX+1));  // k gauss pts ==> accurate for 2k-1 degree polynomial

	Element::PMAX_alpt = AlptBasis::PMAX;	// max polynomial degree for Alpert's basis functions
	Element::PMAX_intp = LagrBasis::PMAX;	// max polynomial degree for interpolation basis functions
	Element::DIM = DIM;			// dimension
	Element::VEC_NUM = VEC_NUM;		// num of unknown variables in PDEs; changed
    Element::is_intp.resize(Element::VEC_NUM);
	for (size_t num = 0; num < Element::VEC_NUM; num++) { Element::is_intp[num] = std::vector<bool>(Element::DIM, true); }

	DGSolution::DIM = Element::DIM;
	DGSolution::VEC_NUM = Element::VEC_NUM;
    DGSolution::ind_var_vec.resize(Element::VEC_NUM);
    for (auto num = 0; num < Element::VEC_NUM; num++) DGSolution::ind_var_vec[num] = num;
    DGAdapt::indicator_var_adapt = DGSolution::ind_var_vec;

	Interpolation::DIM = Element::DIM;
	Interpolation::VEC_NUM = Element::VEC_NUM;
    
	const bool sparse = true;
    const std::string boundary_type = "period";

    // hash key
	Hash hash;

	LagrBasis::set_interp_msh01();
	HermBasis::set_interp_msh01();

    double refine_eps = 1e6;   
	double coarsen_eta = -1;  
	const bool is_adapt_find_ptr_alpt = true;	// variable control if need to adaptively find out pointers related to Alpert basis in DG operators
	const bool is_adapt_find_ptr_intp = true;	// variable control if need to adaptively find out pointers related to interpolation basis in DG operators

    AllBasis<LagrBasis> all_bas_lagr(NMAX);
    AllBasis<HermBasis> all_bas_herm(NMAX);
    AllBasis<AlptBasis> all_bas_alpt(NMAX);

    OperatorMatrix1D<AlptBasis,AlptBasis> oper_matx_alpt(all_bas_alpt, all_bas_alpt, boundary_type);
    OperatorMatrix1D<HermBasis,AlptBasis> oper_matx_herm(all_bas_herm, all_bas_alpt, boundary_type);
    OperatorMatrix1D<LagrBasis,AlptBasis> oper_matx_lagr(all_bas_lagr, all_bas_alpt, boundary_type);

    int N_init = NMAX;



    
    auto init_func_1 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::exp(std::cos(x+2));
        else return std::sin(x*x);
    };
    auto init_func_2 = [](double x, int d)->double 
    {
        assert(d>=0 && d<=1);
        if (d == 0) return std::cos(x+2)/(2.+x);
        else return std::sin(x)*x;
    };
    std::vector<std::function<double(double, int)>> init_func{init_func_1, init_func_2};

    auto func_flux1 = [](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double
    {
        assert(d>=0 && d<=1 && i>=0 && i<=1);
        if (d == 0) 
        {
            if (i == 0)
                return std::sin(u_next[0] * u_next[1])/(1. + std::pow(u_next[0] - 2.*u_next[1], 2.));
            else 
                return std::exp( - u_next[0] * u_next[1]);
        }
        else 
        {
            if (i == 0)
                return std::sin(u_next[0] * u_next[1]);
            else 
                return u_next[0] / (u_next[1] * u_next[1] + 5.);
        }
    };

    auto func_flux2 = [&func_flux1](std::vector<double> u_now, std::vector<double> u_next, int i, int d)->double { return func_flux1(u_next,u_now,i,d); };

    auto func_flux = [&func_flux2](std::vector<double> u, int i, int d)->double { std::vector<double> u_zero(Element::VEC_NUM,0.); return func_flux2(u,u_zero,i,d); };


    // next step u is initialized; current step u is 0
    DGAdapt dg_solu(sparse, N_init, NMAX, all_bas_alpt, all_bas_lagr, all_bas_herm, hash, refine_eps, coarsen_eta, is_adapt_find_ptr_alpt, is_adapt_find_ptr_intp);
    LagrInterpolation interp(dg_solu);    
    dg_solu.set_ucoe_alpt_zero();
    dg_solu.init_separable_system(init_func);
    dg_solu.copy_ucoe_current_to_next();  
    dg_solu.set_ucoe_alpt_zero(); 
    interp.nonlinear_Lagr(func_flux1, Element::is_intp);

    // next step u is 0; current step u is initialized
    DGAdapt dg_solu2(sparse, N_init, NMAX, all_bas_alpt, all_bas_lagr, all_bas_herm, hash, refine_eps, coarsen_eta, is_adapt_find_ptr_alpt, is_adapt_find_ptr_intp);
    LagrInterpolation interp2(dg_solu2);    
    dg_solu2.set_ucoe_alpt_zero();
    dg_solu2.init_separable_system(init_func);
    interp2.nonlinear_Lagr(func_flux2, Element::is_intp);

    double errnorm = 1.;
    Substract_fu_norm(dg_solu,dg_solu2,errnorm);
    EXPECT_NEAR(0.,errnorm,1e-15);

    interp2.nonlinear_Lagr(func_flux, Element::is_intp);
    errnorm = 1.;
    Substract_fu_norm(dg_solu,dg_solu2,errnorm);
    EXPECT_NEAR(0.,errnorm,1e-15);

}


























void Substract_fu_norm(DGSolution & dg_solu1, DGSolution & dg_solu2, double &val)
{
    assert(dg_solu1.dg.size() == dg_solu2.dg.size());
    val = 0;
    for (auto && iter : dg_solu1.dg ) {
        auto iter2 = dg_solu2.dg.find(iter.first);
        assert(iter2 != dg_solu2.dg.end());
        for (auto i = 0; i < Element::VEC_NUM; i++) 
            for (auto j = 0; j < Element::DIM; j++) 
            {
                val += std::fabs( (iter.second.fucoe_intp[i][j] - iter2->second.fucoe_intp[i][j]).norm() );
                val += std::fabs( (iter.second.fp_intp[i][j] - iter2->second.fp_intp[i][j]).norm() );
            }
    }
}






















void Interpolation_New_Test__Lagr_General_Check_Order_System_Separable(int given_PMAX, int given_Lagr_PMAX, int smallNMAX, int largeNMAX, int DIM, int VEC_NUM,
    std::vector<std::function<double(double, int)>> init_func_current, 
    std::vector<std::function<double(double, int)>> init_func_next,
    std::function<double(std::vector<double> /* u_now */, std::vector<double> /* u_next */, int, int)> func_flux, bool thesparse)
{

    // define static variables

    AlptBasis::PMAX = given_PMAX;
    const bool sparse = thesparse;
    LagrBasis::PMAX = given_Lagr_PMAX;

    std::cout << "------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Alpert basis PMAX = " << AlptBasis::PMAX << " , Lagrangian basis PMAX = " << LagrBasis::PMAX << " . Please MANUALLY check the order below." << std::endl;

	LagrBasis::msh_case = 1;

	HermBasis::PMAX = 3;    // Hermite is not tested here though
	HermBasis::msh_case = 1;

    int common_num_gauss_pt = 1+(int)(0.5*(std::max(given_PMAX,given_Lagr_PMAX)+1));  // k gauss pts ==> accurate for 2k-1 degree polynomial

	Element::PMAX_alpt = AlptBasis::PMAX;	// max polynomial degree for Alpert's basis functions
	Element::PMAX_intp = LagrBasis::PMAX;	// max polynomial degree for interpolation basis functions
	Element::DIM = DIM;			// dimension
	Element::VEC_NUM = VEC_NUM;		// num of unknown variables in PDEs; changed
    Element::is_intp.resize(Element::VEC_NUM);
	for (size_t num = 0; num < Element::VEC_NUM; num++) { Element::is_intp[num] = std::vector<bool>(Element::DIM, true); }

	DGSolution::DIM = Element::DIM;
	DGSolution::VEC_NUM = Element::VEC_NUM;
    DGSolution::ind_var_vec.resize(Element::VEC_NUM);
    for (auto num = 0; num < Element::VEC_NUM; num++) DGSolution::ind_var_vec[num] = num;
    DGAdapt::indicator_var_adapt = DGSolution::ind_var_vec;

	Interpolation::DIM = Element::DIM;
	Interpolation::VEC_NUM = Element::VEC_NUM;
    
    const std::string boundary_type = "period";

    // hash key
	Hash hash;

	LagrBasis::set_interp_msh01();
	HermBasis::set_interp_msh01();

    double refine_eps = 1e6;   
	double coarsen_eta = -1;  
	const bool is_adapt_find_ptr_alpt = true;	// variable control if need to adaptively find out pointers related to Alpert basis in DG operators
	const bool is_adapt_find_ptr_intp = true;	// variable control if need to adaptively find out pointers related to interpolation basis in DG operators

    


    std::vector<double> err_fu_prev{0.,0.,0.};
    std::vector<double> err_current_prev{0.,0.,0.};
    std::vector<double> err_next_prev{0.,0.,0.};

    std::vector<std::vector<double>> err_fu2_prev(Element::DIM,{0.,0.,0.});
    std::vector<std::vector<double>> err_fu_DW_prev(Element::DIM,{0.,0.,0.});


    for (int NMAX = smallNMAX; NMAX <= largeNMAX; NMAX=NMAX+1)
    {
        //if (sparse) {
            //std::cout << "When max mesh level N = " << NMAX << ", f(u) order might be reduced to minimum of (interpolation order) " 
            //        << LagrBasis::PMAX+1-(Element::DIM-1.0)*log2(1.0+1.0/NMAX) 
            //        << " and (initial projection order) " << AlptBasis::PMAX+1 << std::endl; //<< AlptBasis::PMAX+1-(Element::DIM-1.0)*log2(NMAX)/NMAX << std::endl;
            // assume error is Delta(h)=|log_2(h)|^{d-1}h^{k+1}, so for h=1/2^N, order = log2(Delta(h)/Delta(h/2)) = k+1-(d-1)*log2(1+1/N)
        //}else {
            std::cout << "When max mesh level N = " << NMAX << ":" << std::endl;
        //}

        AllBasis<LagrBasis> all_bas_lagr(NMAX);
        AllBasis<HermBasis> all_bas_herm(NMAX);
        AllBasis<AlptBasis> all_bas_alpt(NMAX);

        OperatorMatrix1D<AlptBasis,AlptBasis> oper_matx_alpt(all_bas_alpt, all_bas_alpt, boundary_type);
        OperatorMatrix1D<HermBasis,AlptBasis> oper_matx_herm(all_bas_herm, all_bas_alpt, boundary_type);
        OperatorMatrix1D<LagrBasis,AlptBasis> oper_matx_lagr(all_bas_lagr, all_bas_alpt, boundary_type);

        int N_init = NMAX;
        DGAdapt dg_solu(sparse, N_init, NMAX, all_bas_alpt, all_bas_lagr, all_bas_herm, hash, refine_eps, coarsen_eta, is_adapt_find_ptr_alpt, is_adapt_find_ptr_intp);

        LagrInterpolation interp(dg_solu);


        // Initial condition, projection only, for next step u        
        dg_solu.set_ucoe_alpt_zero();
        dg_solu.init_separable_system(init_func_next);
        const int num_gauss_pt_1 = common_num_gauss_pt;
        std::vector<double> err_next = dg_solu.get_error_separable_system(init_func_next, num_gauss_pt_1);
        std::cout << "L1, L2 and Linf error for next step u projection: " << err_next[0] << ", " << err_next[1] << ", " << err_next[2];
        if (NMAX > smallNMAX) {
            std::cout << ", orders are " << log2(err_next_prev[0]/err_next[0]) << ", " << log2(err_next_prev[1]/err_next[1]) << ", " << log2(err_next_prev[2]/err_next[2]);
        }
        std::cout << std::endl;
        dg_solu.copy_ucoe_current_to_next();



        
        // Initial condition, projection only, for current u
        ///////////////////////// Imporant notice: when init_separable_system is called (which called init_elem_separable inside) ////////////////////////////
        ///////////////////////// the code !!!!!!!!!!! ADD !!!!!!!!!!!!! the given initial condition to current ucoe_alpt ////////////////////////////////////
        ///////////////////////// so if you do not want to add, be sure to clean up ucoe_alpt first //////////////////////////////////////////////////////////
        dg_solu.set_ucoe_alpt_zero();
        dg_solu.init_separable_system(init_func_current);
        const int num_gauss_pt_2 = common_num_gauss_pt;
        std::vector<double> err_current = dg_solu.get_error_separable_system(init_func_current, num_gauss_pt_2);
        std::cout << "L1, L2 and Linf error for current step u projection: " << err_current[0] << ", " << err_current[1] << ", " << err_current[2];
        if (NMAX > smallNMAX) {
            std::cout << ", orders are " << log2(err_current_prev[0]/err_current[0]) << ", " << log2(err_current_prev[1]/err_current[1]) << ", " << log2(err_current_prev[2]/err_current[2]);
        }
        std::cout << std::endl;




        // Flux function f(u_current,u_next) Lagrangian interpolation
        interp.nonlinear_Lagr(func_flux, Element::is_intp);



        // test accuracy of f(u_now,u_next) interpolation
        // init_func_1 for next step, init_func_2 for current step
        auto fu_ext = [&init_func_current,&init_func_next,&func_flux](std::vector<double> x, int i, int d)->double 
        {
            assert(i>=0 && i<Element::VEC_NUM && d>=0 && d<Element::DIM && x.size()==Element::DIM);
            //std::vector<double> u_now = {init_func_current[0](x[0],0), init_func_current[1](x[0],0), init_func_current[2](x[0],0)};
            //std::vector<double> u_next = {init_func_next[0](x[0],0), init_func_next[1](x[0],0), init_func_next[2](x[0],0)};
            std::vector<double> u_now(Element::VEC_NUM,1.);
            std::vector<double> u_next(Element::VEC_NUM,1.);
            for (auto li = 0; li < Element::DIM; li++) {
                for (auto lj = 0; lj < Element::VEC_NUM; lj++) {
                    u_now[lj] *= init_func_current[lj](x[li],li);
                    u_next[lj] *= init_func_next[lj](x[li],li);
                }
            }
            return func_flux(u_now,u_next,i,d);
        };

        const int num_gauss_pt = common_num_gauss_pt;
        std::vector<double> err_fu = dg_solu.get_error_Lag(fu_ext, num_gauss_pt, Element::is_intp);
        std::cout << "L1, L2 and Linf error for f(u): " << err_fu[0] << ", " << err_fu[1] << ", " << err_fu[2];
        if (NMAX > smallNMAX) {
            std::cout << ", orders are " << log2(err_fu_prev[0]/err_fu[0]) << ", " << log2(err_fu_prev[1]/err_fu[1]) << ", " << log2(err_fu_prev[2]/err_fu[2]);
        }
        std::cout << std::endl;


        // test accuracy of f(u) by dimension using info in fp_intp (info at interpolation points)
        //std::vector<std::vector<std::vector<double>>> err_fu_CDW = dg_solu.get_error_Lag_CompDimWise(fu_ext, num_gauss_pt, Element::is_intp);
        std::vector<std::vector<double>> err_fu_DW = dg_solu.get_error_Lag_DimWise(fu_ext, num_gauss_pt, Element::is_intp);
        for (auto dd = 0; dd < Element::DIM; dd++)
        {
            // for (auto ii = 0; ii < Element::VEC_NUM; ii++)
            //     if (Element::is_intp[ii][dd])
            //     {
            //         std::cout << "L1, L2, Linf error for f(u) on component " << ii << " and dimension " << dd << " are " 
            //                   << err_fu_CDW[ii][dd][0] << ", " << err_fu_CDW[ii][dd][1] << ", " << err_fu_CDW[ii][dd][2] << std::endl;
            //     }
            bool is_intp_DIM = false;
            for (auto ii = 0; ii < Element::VEC_NUM; ii++) is_intp_DIM = (is_intp_DIM || Element::is_intp[ii][dd]);
            if (is_intp_DIM) {
                std::cout << "L1, L2, Linf error for f(u) on dimension " << dd << ": " 
                          << err_fu_DW[dd][0] << ", " << err_fu_DW[dd][1] << ", " << err_fu_DW[dd][2];
                if (NMAX>smallNMAX) {
                    std::cout << ", orders are " << log2(err_fu_DW_prev[dd][0]/err_fu_DW[dd][0]) << ", " << log2(err_fu_DW_prev[dd][1]/err_fu_DW[dd][1]) << ", " << log2(err_fu_DW_prev[dd][2]/err_fu_DW[dd][2]);
                }
                std::cout << std::endl;
                err_fu_DW_prev[dd] = err_fu_DW[dd];
            }
        }


        // test accuracy of f(u) by dimension using info in fucoe_intp (alpert basis coefficients)
        for (auto dd = 0; dd < Element::DIM; dd++) 
        {
            std::vector<std::function<double(std::vector<double>)>> fu_ext_dd(Element::VEC_NUM);
            for (auto ii = 0; ii < Element::VEC_NUM; ii++)
                fu_ext_dd[ii] = [ii,dd,&fu_ext](std::vector<double> x)->double { return fu_ext(x,ii,dd); }; // take care of the capture here, especially ii
            // size of iter.second.ucoe_alpt[ii] is (AlptBasis::PMAX+1)^2, while size of iter.second.fucoe_intp[ii][dd] is (LagrBasis::PMAX+1)^2
            // Element::is_intp should be all true; DGSolution::ind_var_vec and DGAdapt::indicator_var_adapt should be full; otherwise there might be problem here
            BilinearFormIntp interp_bf(dg_solu, oper_matx_lagr, oper_matx_herm);
            for (auto ii = 0; ii < Element::VEC_NUM; ii++)
                interp_bf.assemble_matrix_intp(1, dd, interp_bf.oper_matx_lagr_ptr->u_v, interp_bf.oper_matx_lagr_ptr->u_v, "vol", ii, ii);
            ODESolver odesolver(dg_solu);
            odesolver.fucoe_to_eigenvec(dd);
            odesolver.ucoe = interp_bf.mat * odesolver.fucoe;
            odesolver.eigenvec_to_ucoe();
            
            std::vector<double> err_fu2 = dg_solu.get_error_no_separable_system(fu_ext_dd,common_num_gauss_pt);
            std::cout << "L1, L2 and Linf error for f(u) (alpert basis) on dimension " << dd << ": " << err_fu2[0] << ", " << err_fu2[1] << ", " << err_fu2[2];
            if (NMAX > smallNMAX) {
                std::cout << ", orders are " << log2(err_fu2_prev[dd][0]/err_fu2[0]) << ", " << log2(err_fu2_prev[dd][1]/err_fu2[1]) << ", " << log2(err_fu2_prev[dd][2]/err_fu2[2]);
            }
            std::cout << std::endl;
            err_fu2_prev[dd] = err_fu2;
            
            if (NMAX == smallNMAX) {
                auto fu_ext_pt2norm = [&fu_ext_dd](std::vector<double> x)->double { 
                    double ret=0.; 
                    for (auto ii=0;ii<Element::VEC_NUM;ii++) 
                        ret += mypow(fu_ext_dd[ii](x),2);
                    return std::sqrt(ret);
                };
                Quad fu_ext_norm_cal(Element::DIM);
                std::vector<double> fu_ext_norm = fu_ext_norm_cal.norm_multiD(fu_ext_pt2norm,NMAX,AlptBasis::PMAX+1);
                std::cout << "At dimension " << dd << " flux norm is " << fu_ext_norm[0] << ", " << fu_ext_norm[1] << ", " << fu_ext_norm[2] << std::endl;
            }
        }


        err_current_prev = err_current;
        err_next_prev = err_next;
        err_fu_prev = err_fu;
    }

    std::cout << "------------------------------------------------------------------------------------" << std::endl;
    
}









































