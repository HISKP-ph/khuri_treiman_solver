#include "array.h"
#include "constants.h"
#include<fstream>

using namespace array;
using constants::pi;




Array::Array(double mass_1, double mass_2, double mass_3, double mass_4, double cutoff)
:
s_I{pow(mass_1+mass_2,2)}, s_II{mass_1*(std::pow(mass_4,2.)-std::pow(mass_3,2.))/(mass_3+mass_2)},
s_III{pow(mass_3-mass_4,2)}, s_IV{pow(mass_4+mass_3,2)},
t_I{std::pow(mass_1+mass_3,2)},
t_III{std::pow(mass_4-mass_2,2)}, t_IV{std::pow(mass_4+mass_2,2)},
decay_mass{mass_4},
cutoff{cutoff}
{}

ArrayEta3Pi::ArrayEta3Pi(double mass_1, double mass_4, double cutoff)
:
Array{mass_1,mass_1,mass_1,mass_4,cutoff}
{create_interval();}

ArrayEtap3Pi::ArrayEtap3Pi(double mass_1, double mass_4, double cutoff)
:
Array{mass_1,mass_1,mass_1,mass_4,cutoff}
{create_interval();}


ArrayEtapEtaPiPi::ArrayEtapEtaPiPi(double mass_1, double mass_3, double mass_4, double cutoff)
:
Array{mass_1,mass_1,mass_3,mass_4,cutoff}
{create_interval();}

ArrayV3Pi::ArrayV3Pi(double mass_1, double mass_4, double cutoff)
:
Array{mass_1,mass_1,mass_1,mass_4,cutoff}
{create_interval();}

ArrayX3Pi::ArrayX3Pi(double mass_1, double mass_4, double cutoff)
:
Array{mass_1,mass_1,mass_1,mass_4,cutoff}
{create_interval();}

ArrayT3Pi::ArrayT3Pi(double mass_1, double mass_4, double cutoff)
:
Array{mass_1,mass_1,mass_1,mass_4,cutoff}
{create_interval();}


void Array::asymptotic_spacing(std::vector<double>& list, double min, double max){
    double value=min;
    // narrow spacing slightly above min
    for (double j=eps; value<=min+0.1 ; j*=2){
        value+=j;
        list.push_back(value);
    }
    
    // constant spacing between min and max
    for (value=min+0.2; value<max-0.3 ; value+=0.1){
        list.push_back(value);
    }
    if(max-0.3 > value){
        list.push_back(max-0.3);
    }
    if(max-0.27 > value){
        list.push_back(max-0.27);
    }
    value=max-0.2;


    
    // narrow spacing slightly below max
    for (double j=0.1; value<max-eps ; j*=0.5){
        value+=j;
        list.push_back(value);
    }
    return;
}


void Array::constant_spacing(std::vector<double>& list, double min, double max, double step_size){
    // constant spacing between min and max
    for (double value=min; value<=max; value+=step_size){
        list.push_back(value);
    }
    return;
}

void ArrayEta3Pi::create_interval(){
    constant_spacing(interval_s, -cutoff,-1,2);
    constant_spacing(interval_disp, -cutoff,-1,2);
    
    interval_s.push_back(-eps);
    interval_disp.push_back(-eps);
    
    interval_s.push_back(eps);
    interval_disp.push_back(eps);
    
    
    asymptotic_spacing(interval_s, .1-eps, s_I);
    asymptotic_spacing(interval_disp, .1-eps, s_I);
    
    
    asymptotic_spacing(interval_s, s_I, s_II);
    asymptotic_spacing(interval_s_th, s_I, s_II);
    asymptotic_spacing(interval_disp, s_I, s_II);
    
    asymptotic_spacing(interval_s, s_II, s_III);
    asymptotic_spacing(interval_s_th, s_II, s_III);
   
    asymptotic_spacing(interval_s, s_III, s_IV);
    asymptotic_spacing(interval_s_th, s_III, s_IV);
    asymptotic_spacing(interval_egg, s_III, s_IV);

    
    double large=std::pow(decay_mass+2,2);
    if (large < 80)
    {
        large = 80;
    }

    constant_spacing(interval_s, s_IV+0.2, large, .2+eps);
    constant_spacing(interval_s_th, s_IV+0.2, large, .2+eps);


    constant_spacing(interval_s, large, cutoff, 2);
    constant_spacing(interval_s_th, large,cutoff, 2);

    double pseudo_thr;
    if (decay_mass < 15){
        pseudo_thr = 7;
    }
    else{
        pseudo_thr = decay_mass+3.;
    }
    constant_spacing(interval_disp, s_II+0.2, s_III-pseudo_thr, 0.1);
    constant_spacing(interval_disp, s_III+pseudo_thr, s_IV-2.5, 0.1);
    constant_spacing(interval_disp, s_IV+2.5, large, 0.1);
    constant_spacing(interval_disp, large, cutoff-1e-6, 2);

    constant_spacing(interval_disp_expansion, s_III-2., s_III+2., 0.01);

    for (double angle=0.0000000000001; angle<= 2.*pi(); angle+=0.01){
        interval_phi.push_back(angle);
    }
    interval_phi.push_back(2.*pi()-0.0000000000001);
    return;
}

// old eta3pi not used in current project
// void ArrayEta3Pi::create_interval(){
//     double cusp=4.0*pow(constants::mass_kaon(),2.0)/pow(constants::mass_pi(),2.0);
//     constant_spacing(interval_s, -cutoff,-1,2);
//     constant_spacing(interval_disp, -cutoff,-1,2);
    
//     interval_s.push_back(-eps);
//     interval_disp.push_back(-eps);
    
//     interval_s.push_back(eps);
//     interval_disp.push_back(eps);
    
    
//     asymptotic_spacing(interval_s, .1-eps, s_I);
//     asymptotic_spacing(interval_disp, .1-eps, s_I);
    
    
//     asymptotic_spacing(interval_s, s_I, s_II);
//     asymptotic_spacing(interval_s_th, s_I, s_II);
//     asymptotic_spacing(interval_disp, s_I, s_II);
    
//     asymptotic_spacing(interval_s, s_II, s_III);
//     asymptotic_spacing(interval_s_th, s_II, s_III);
   
//     asymptotic_spacing(interval_s, s_III, s_IV);
//     asymptotic_spacing(interval_s_th, s_III, s_IV);
//     asymptotic_spacing(interval_egg, s_III, s_IV);

//     double large=std::pow(decay_mass+2,2);
//     if (large < 80)
//     {
//         large = 80;
//     }

//     constant_spacing(interval_s, s_IV+0.2, large, .2+eps);
//     constant_spacing(interval_s_th, s_IV+0.2, large, .2+eps);


//     constant_spacing(interval_s, large, cutoff, 2);
//     constant_spacing(interval_s_th, large,cutoff, 2);

//     //special for each decay
//     constant_spacing(interval_disp, s_II+0.2, s_III-.9, 0.1);
//     constant_spacing(interval_disp, s_III+.9, s_IV-2.5, 0.1);
//     constant_spacing(interval_disp, s_IV+2.5, cusp-.1, 0.1);
//     interval_disp.push_back(cusp);
//     constant_spacing(interval_disp, cusp+.1, large, 0.1);
//     // constant_spacing(interval_disp, s_IV+2.5, large, 0.1); //needed to calculate for decay_mass -> large, then only P-waves ok
//     constant_spacing(interval_disp, large, cutoff-1e-6, 2);

//     // Lists used to expand the interval. Only evaluate the expansion along these short lists to save computation time.
//     constant_spacing(interval_disp_expansion, s_III-2., s_III+2., 0.01);

//     for (double angle=0.0000000000001; angle<= 2.*pi(); angle+=0.01){
//         interval_phi.push_back(angle);
//     }
//     interval_phi.push_back(2.*pi()-0.0000000000001);
//     return;
// }

void ArrayEtap3Pi::create_interval(){
    double cusp=4.0*pow(constants::mass_kaon(),2.0)/pow(constants::mass_pi(),2.0);
    constant_spacing(interval_s, -cutoff,-1,2);
    constant_spacing(interval_disp, -cutoff,-1,2);
    
    interval_s.push_back(-eps);
    interval_disp.push_back(-eps);
    
    interval_s.push_back(eps);
    interval_disp.push_back(eps);
    
    
    asymptotic_spacing(interval_s, .1-eps, s_I);
    asymptotic_spacing(interval_disp, .1-eps, s_I);
    
    
    asymptotic_spacing(interval_s, s_I, s_II);
    asymptotic_spacing(interval_s_th, s_I, s_II);
    asymptotic_spacing(interval_disp, s_I, s_II);
    
    asymptotic_spacing(interval_s, s_II, s_III);
    asymptotic_spacing(interval_s_th, s_II, s_III);
   
    asymptotic_spacing(interval_s, s_III, s_IV);
    asymptotic_spacing(interval_s_th, s_III, s_IV);
    asymptotic_spacing(interval_egg, s_III, s_IV);

    double large=std::pow(decay_mass+1.5,2);

    constant_spacing(interval_s, s_IV+0.2, large, .2+eps);
    constant_spacing(interval_s_th, s_IV+0.2, large, .2+eps);


    constant_spacing(interval_s, large, cutoff, 2);
    constant_spacing(interval_s_th, large,cutoff, 2);

    //special for each decay
    constant_spacing(interval_disp, s_II+0.2, s_III-7, 0.1);
    constant_spacing(interval_disp, s_III+7, cusp-.5, 0.1);
    interval_disp.push_back(cusp);
    constant_spacing(interval_disp, cusp+.1,  s_IV-.3, 0.1);
    constant_spacing(interval_disp, s_IV+.5, large, 0.1);
    constant_spacing(interval_disp, large, cutoff-1e-6, 2);
    
    // Lists used to expand the interval. Only evaluate the expansion along these short lists to save computation time.
    constant_spacing(interval_disp_expansion, s_III-2., s_III+2., 0.01);

    for (double angle=0.0000000000001; angle<= 2.*pi(); angle+=0.01){
        interval_phi.push_back(angle);
    }
    interval_phi.push_back(2.*pi()-0.0000000000001);
    return;
}


void ArrayEtapEtaPiPi::create_interval(){
    double cusp=4.0*pow(constants::mass_kaon(),2.0)/pow(constants::mass_pi(),2.0);
    constant_spacing(interval_s, -cutoff,-1,2);
    constant_spacing(interval_disp, -cutoff,-1,2);
    
    interval_s.push_back(-eps);
    interval_disp.push_back(-eps);
    
    interval_s.push_back(eps);
    interval_disp.push_back(eps);
    
    
    asymptotic_spacing(interval_s, .1-eps, s_I);
    asymptotic_spacing(interval_disp, .1-eps, s_I);
    
    
    asymptotic_spacing(interval_s, s_I, s_II);
    asymptotic_spacing(interval_s_th, s_I, s_II);
    asymptotic_spacing(interval_disp, s_I, s_II);
    
    asymptotic_spacing(interval_s, s_II, s_III);
    asymptotic_spacing(interval_s_th, s_II, s_III);
   
    asymptotic_spacing(interval_s, s_III, s_IV);
    asymptotic_spacing(interval_s_th, s_III, s_IV);

    
    // use a narrow spacing up to some large value, from there on use less points
    double large=std::pow(decay_mass+1.5,2);

    constant_spacing(interval_s, s_IV+0.2, large, .2+eps);
    constant_spacing(interval_s_th, s_IV+0.2, large, .2+eps);


    constant_spacing(interval_s, large, cutoff, 2);
    constant_spacing(interval_s_th, large,cutoff, 2);

    //special for each decay
    constant_spacing(interval_disp, s_II+0.2, s_III-2, 0.1);
    constant_spacing(interval_disp, s_III+2, cusp-.1, 0.1);
    interval_disp.push_back(cusp);
    constant_spacing(interval_disp, cusp+.1,  s_IV-2, 0.1);
    constant_spacing(interval_disp, s_IV+2, large, 0.1);
    constant_spacing(interval_disp, large, cutoff-1e-6, 2);

    // Lists used to expand the interval. Only evaluate the expansion along these short lists to save computation time.
    constant_spacing(interval_disp_expansion, s_III-2., s_III+2., 0.01);
    constant_spacing(interval_t_disp_expansion, t_III-2., t_III+2., 0.01);




    // here create the list that is only used for eta'-> eta pi pi

    // All values for the curve parameters should be included in this small array
    constant_spacing(interval_y,0.7,3.3,0.05);


    constant_spacing(interval_t, -cutoff,-1,2);
    interval_t.push_back(-eps);
    interval_t.push_back(eps);
    asymptotic_spacing(interval_t, .1-eps, t_I-.1);

    asymptotic_spacing(interval_t, t_I+.1, t_III);
    asymptotic_spacing(interval_t_th, t_I, t_III);

    asymptotic_spacing(interval_t, t_III, t_IV);
    asymptotic_spacing(interval_t_th, t_III, t_IV);
    asymptotic_spacing(interval_t_path_spline, t_III, t_IV);

    constant_spacing(interval_t, t_IV+0.2, large, .2+eps);
    constant_spacing(interval_t_th, t_IV+0.2, large, .2+eps);

    constant_spacing(interval_t, large, cutoff, 2);
    constant_spacing(interval_t_th, large,cutoff, 2);



    //-- special treatment of list to evaluate dispersionintegral ------------
    constant_spacing(interval_t_disp, -cutoff,-1,2);
    interval_t_disp.push_back(-eps);
    interval_t_disp.push_back(eps);   
    constant_spacing(interval_t_disp, .1+eps, s_III-1,.1); 
    constant_spacing(interval_t_disp, s_III+1, t_I,.1); 
    constant_spacing(interval_t_disp, t_I+0.1, t_III-3.5, 0.1);
    constant_spacing(interval_t_disp, t_III+3.5, cusp-.1, 0.1);
    interval_t_disp.push_back(cusp);
    constant_spacing(interval_t_disp, cusp+.1,  t_IV-1.5, 0.1);
    constant_spacing(interval_t_disp, t_IV+1.5, large, 0.1);
    constant_spacing(interval_t_disp, large, cutoff-1e-6, 2);
    
    return;
}

void ArrayV3Pi::create_interval(){
    constant_spacing(interval_s, -cutoff,-1,2);
    constant_spacing(interval_disp, -cutoff,-1,2);
    
    interval_s.push_back(-eps);
    interval_disp.push_back(-eps);
    
    interval_s.push_back(eps);
    interval_disp.push_back(eps);
    
    
    asymptotic_spacing(interval_s, .1-eps, s_I);
    asymptotic_spacing(interval_disp, .1-eps, s_I);
    
    
    asymptotic_spacing(interval_s, s_I, s_II);
    asymptotic_spacing(interval_s_th, s_I, s_II);
    asymptotic_spacing(interval_disp, s_I, s_II);
    
    asymptotic_spacing(interval_s, s_II, s_III);
    asymptotic_spacing(interval_s_th, s_II, s_III);
   
    asymptotic_spacing(interval_s, s_III, s_IV);
    asymptotic_spacing(interval_s_th, s_III, s_IV);
    asymptotic_spacing(interval_egg, s_III, s_IV);

    
    double large=std::pow(decay_mass+2,2);
    if (large < 80)
    {
        large = 80;
    }

    constant_spacing(interval_s, s_IV+0.2, large, .2+eps);
    constant_spacing(interval_s_th, s_IV+0.2, large, .2+eps);


    constant_spacing(interval_s, large, cutoff, 2);
    constant_spacing(interval_s_th, large,cutoff, 2);

    //special for each decay
    constant_spacing(interval_disp, s_II+0.2, s_III-7, 0.1);
    constant_spacing(interval_disp, s_III+7, s_IV-2.5, 0.1);
    constant_spacing(interval_disp, s_IV+2.5, large, 0.1);
    constant_spacing(interval_disp, large, cutoff-1e-6, 2);

    // Lists used to expand the interval. Only evaluate the expansion along these short lists to save computation time.
    constant_spacing(interval_disp_expansion, s_III-2., s_III+2., 0.01);

    for (double angle=0.0000000000001; angle<= 2.*pi(); angle+=0.01){
        interval_phi.push_back(angle);
    }
    interval_phi.push_back(2.*pi()-0.0000000000001);
    return;
}

void ArrayX3Pi::create_interval(){
    constant_spacing(interval_s, -cutoff,-1,2);
    constant_spacing(interval_disp, -cutoff,-1,2);
    
    interval_s.push_back(-eps);
    interval_disp.push_back(-eps);
    
    interval_s.push_back(eps);
    interval_disp.push_back(eps);
    
    
    asymptotic_spacing(interval_s, .1-eps, s_I);
    asymptotic_spacing(interval_disp, .1-eps, s_I);
    
    
    asymptotic_spacing(interval_s, s_I, s_II);
    asymptotic_spacing(interval_s_th, s_I, s_II);
    asymptotic_spacing(interval_disp, s_I, s_II);
    
    asymptotic_spacing(interval_s, s_II, s_III);
    asymptotic_spacing(interval_s_th, s_II, s_III);
   
    asymptotic_spacing(interval_s, s_III, s_IV);
    asymptotic_spacing(interval_s_th, s_III, s_IV);
    asymptotic_spacing(interval_egg, s_III, s_IV);

    
    double large=std::pow(decay_mass+2,2);
    if (large < 80)
    {
        large = 80;
    }

    constant_spacing(interval_s, s_IV+0.2, large, .2+eps);
    constant_spacing(interval_s_th, s_IV+0.2, large, .2+eps);


    constant_spacing(interval_s, large, cutoff, 2);
    constant_spacing(interval_s_th, large,cutoff, 2);

    //special for each decay
    constant_spacing(interval_disp, s_II+0.2, s_III-.9, 0.1);
    constant_spacing(interval_disp, s_III+.9, s_IV-2.5, 0.1);
    constant_spacing(interval_disp, s_IV+2.5, large, 0.1);
    constant_spacing(interval_disp, large, cutoff-1e-6, 2);

    // Lists used to expand the interval. Only evaluate the expansion along these short lists to save computation time.
    constant_spacing(interval_disp_expansion, s_III-2., s_III+2., 0.01);

    for (double angle=0.0000000000001; angle<= 2.*pi(); angle+=0.01){
        interval_phi.push_back(angle);
    }
    interval_phi.push_back(2.*pi()-0.0000000000001);
    return;
}

void ArrayT3Pi::create_interval(){
    constant_spacing(interval_s, -cutoff,-1,2);
    constant_spacing(interval_disp, -cutoff,-1,2);
    
    interval_s.push_back(-eps);
    interval_disp.push_back(-eps);
    
    interval_s.push_back(eps);
    interval_disp.push_back(eps);
    
    
    asymptotic_spacing(interval_s, .1-eps, s_I);
    asymptotic_spacing(interval_disp, .1-eps, s_I);
    
    
    asymptotic_spacing(interval_s, s_I, s_II);
    asymptotic_spacing(interval_s_th, s_I, s_II);
    asymptotic_spacing(interval_disp, s_I, s_II);
    
    asymptotic_spacing(interval_s, s_II, s_III);
    asymptotic_spacing(interval_s_th, s_II, s_III);
   
    asymptotic_spacing(interval_s, s_III, s_IV);
    asymptotic_spacing(interval_s_th, s_III, s_IV);
    asymptotic_spacing(interval_egg, s_III, s_IV);

    
    double large=std::pow(decay_mass+2,2);
    if (large < 80)
    {
        large = 80;
    }

    constant_spacing(interval_s, s_IV+0.2, large, .2+eps);
    constant_spacing(interval_s_th, s_IV+0.2, large, .2+eps);


    constant_spacing(interval_s, large, cutoff, 2);
    constant_spacing(interval_s_th, large,cutoff, 2);

    //special for each decay
    if (decay_mass < 8)
    {
        constant_spacing(interval_disp, s_II+0.2, s_III-7, 0.1);
        constant_spacing(interval_disp, s_III+7, s_IV-2.5, 0.1);
    }
    else if (decay_mass < 14){
        constant_spacing(interval_disp, s_II+0.2, s_III-10, 0.1);
        constant_spacing(interval_disp, s_III+10, s_IV-2.5, 0.1);   
    }
    else{
        constant_spacing(interval_disp, s_II+0.2, s_III-15, 0.1);
        constant_spacing(interval_disp, s_III+15, s_IV-2.5, 0.1);   
    }
    constant_spacing(interval_disp, s_IV+2.5, large, 0.1);
    constant_spacing(interval_disp, large, cutoff-1e-6, 2);

    // Lists used to expand the interval. Only evaluate the expansion along these short lists to save computation time.
    constant_spacing(interval_disp_expansion, s_III-2., s_III+2., 0.01);

    for (double angle=0.0000000000001; angle<= 2.*pi(); angle+=0.01){
        interval_phi.push_back(angle);
    }
    interval_phi.push_back(2.*pi()-0.0000000000001);
    return;
}
