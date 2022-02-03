//
//  RunIntervention.h
//  
//
//  Created by David Hodgson on 05/07/2021.
//

#ifndef RunIntervention_h
#define RunIntervention_h

#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include "ascent/Ascent.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]
using namespace Rcpp;
using namespace std;
using namespace Eigen;
using namespace asc;

// Stuff for random number generation
std::random_device dev;
std::mt19937 engine(dev());
typedef boost::mt19937 PRNG_s;
PRNG_s rng(engine());

class RunInterventions {
public:
    int A;

    vector<double > p_mat;
    vector<double > u_p;

    MatrixXd cnt_matrix_p;
    MatrixXd cnt_matrix_p_h;
    MatrixXd pwp_p;
    MatrixXd pwn_p;
    MatrixXd nwp_p;
    MatrixXd nwn_p;
    
    MatrixXd cnt_matrix_c;
    MatrixXd cnt_matrix_c_h;
    MatrixXd pwp_c;
    MatrixXd pwn_c;
    MatrixXd nwp_c;
    MatrixXd nwn_c;
    
    List SampleWeekly(List vac_calendar, List vac_dose, double cov_c, List vac_info, VectorXd posteriors);
    List SampleMonthly(List vac_calendar, List vac_dose, double cov_c, List vac_info, VectorXd posteriors);
    NumericMatrix getDebug(List vac_calendar, List vac_dose, double cov_c, List vac_info, VectorXd currentParamValues);

    /**********************************************************************/
    /** CwX Contacts for persons who are cocooned ****/
    /**********************************************************************/
    
    vector<vector<double > > get_cwn(double prop_c, char s)
    {
        vector<vector<double > > cwn_e(A,vector<double>(A, 0));
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (nwn_p(i, j));
            
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwn_e[i][j] = (nwn_p(i, j))*(1-prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (pwn_p(i, j));
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwn_e[i][j] = (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(1-prop_c);

        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (nwn_c(i, j));
            
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwn_e[i][j] = (nwn_c(i, j))*(1-prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (pwn_c(i, j));
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwn_e[i][j] = (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return cwn_e;
    }
    
    vector<vector<double > > get_cwp(double prop_c, char s)
    {
        vector<vector<double > > cwp_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwp_e[i][j] = (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(p_mat[j])*(1-prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwp_e[i][j] = (pwp_p(i, j))*(1-prop_c); // parental w. infant plus non-parental with infant who are cocooned
    
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwp_e[i][j] = (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(p_mat[j])*(1-prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwp_e[i][j] = (pwp_c(i, j))*(1-prop_c); // parental w. infant plus non-parental with infant who are cocooned
        }
        else
            cout << "Error cont" << endl;
        
        return cwp_e;
    }
    
    vector<vector<double > > get_cwc(double prop_c, char s)
    {
        vector<vector<double > > cwc_e(A,vector<double>(A,0));

        if (s == 'p')
        {
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwc_e[i][j] = (nwn_p(i, j))*(prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwc_e[i][j] = cnt_matrix_p_h(i, j) + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(prop_c);
            
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwc_e[i][j] = cnt_matrix_p_h(i, j)/2.0 + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(p_mat[j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwc_e[i][j] = (pwp_p(i, j))*(prop_c); // parental w. infant plus non-parental with infant who are cocooned
    
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwc_e[i][j] = (nwn_c(i, j))*(prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwc_e[i][j] = cnt_matrix_c_h(i, j) + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(prop_c);
            
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwc_e[i][j] = cnt_matrix_c_h(i, j)/2.0 + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(p_mat[j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwc_e[i][j] = (pwp_c(i, j))*(prop_c); // parental w. infant plus non-parental with infant who are cocooned

        }
        else
            cout << "Error cont" << endl;
        return cwc_e;
    }
    
    /**********************************************************************/
    /** PwX Contacts for persons who are mothers but not cocooned ****/
    /**********************************************************************/
    
    
    vector<vector<double > > get_pwn(double prop_c, char s)
    {
        vector<vector<double > > pwn_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < A; j++) // all
                    pwn_e[i][j] = pwn_p(i, j); // all contacts
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwn_e[i][j] = cnt_matrix_p_h(i, j) + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(1-prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < A; j++) // all
                    pwn_e[i][j] = pwn_c(i, j); // all contacts
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwn_e[i][j] = cnt_matrix_c_h(i, j) + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return pwn_e;
    }
    
    vector<vector<double > > get_pwp(double prop_c, char s)
    {
        vector<vector<double > > pwp_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // all
                    pwp_e[i][j] = pwp_p(i, j)*(1-prop_c); // all contacts
        }
        else if (s == 'c')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // all
                    pwp_e[i][j] = pwp_c(i, j)*(1-prop_c); // all contacts
        }
        else
            cout << "Error cont" << endl;
        
        return pwp_e;
    }
    
    vector<vector<double > > get_pwc(double prop_c, char s)
    {
        // Non cocooned parent -> cocooned infant (family and outside family)
        vector<vector<double > > pwc_e(A,vector<double>(A,0));
        
        if (s == 'p')
        { 
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwc_e[i][j] = (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    pwc_e[i][j] = pwn_p(i, j)*prop_c; //
        }
        else if (s == 'c')
        {
            // Parent -> parents (all)
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwc_e[i][j] = (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    pwc_e[i][j] = pwn_c(i, j)*prop_c; //
        }
        else
            cout << "Error cont" << endl;
        
        return pwc_e;
    }
    
    
 /**********************************************************************/
/** NwX Contacts for persons who are neither cocooned nor mothers ****/
 /**********************************************************************/
    
    vector<vector<double > > get_nwn(double prop_c, char s)
    {
        // cnt_matrix_p cnt_matrix_p_h
        vector<vector<double > > nwn_e(A,vector<double>(A,0));
        if (s == 'p')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < A; j++) //     infants
                    nwn_e[i][j] = nwn_p(i, j); //number of contacts made by population with non-cocooned prop_ci
            
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwn_e[i][j] = nwn_p(i, j)*(1-prop_c); //number
        }
        else if (s == 'c')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < A; j++) //     infants
                    nwn_e[i][j] = nwn_c(i, j); //number of contacts made by population with non-cocooned prop_ci
            
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwn_e[i][j] = nwn_c(i, j)*(1-prop_c); //numbe

        }
        else
            cout << "Error cont" << endl;

        return nwn_e;
    }
    
    vector<vector<double > > get_nwp(double prop_c, char s)
    {
        // Population with non-cocooned parents)
        vector<vector<double > > nwp_e(A,vector<double>(A,0));
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = (cnt_matrix_p_h(i, j)/2.0) + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*p_mat[j]*(1-prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = cnt_matrix_p(i, j)*p_mat[j]*(1-prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = (cnt_matrix_c_h(i, j)/2.0) + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*p_mat[j]*(1-prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = cnt_matrix_c(i, j)*p_mat[j]*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return nwp_e;
    }
    
    vector<vector<double > > get_nwc(double prop_c, char s)
    {
        
        vector<vector<double > > nwc_e(A,vector<double>(A,0));
        if (s == 'p')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwc_e[i][j] = nwn_p(i, j)*(prop_c); //numbe
            
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*p_mat[j]*(prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = cnt_matrix_p(i, j)*p_mat[j]*(prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwc_e[i][j] = nwn_c(i, j)*(prop_c); //numbe
            
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*p_mat[j]*(prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = cnt_matrix_c(i, j)*p_mat[j]*(prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return nwc_e;
    }
    
    vector< double >  populationPerAgeGroup, eta, modelIncidencePerTime, pA, ep_t;
    double dailyBirthRate, totPopulation;
    NumericVector ageStratification;
    
    int dayNoAfterBurn, weekNo, monthNo;
    double valueLogLikelihood;
    
    // Defined later but of size A
    double currentODETime;
    double run_start, run_burn, run_full, dt;

    RunInterventions(double dailyBirthRate_t, double totPopulation_t, NumericVector ageStratification_t): dailyBirthRate(dailyBirthRate_t), totPopulation(totPopulation_t), ageStratification(ageStratification_t){
           
            A = ageStratification.size();
            eta.push_back(0);
            for (int i = 0; i < A-1; i++){
                populationPerAgeGroup.push_back(dailyBirthRate*365*(ageStratification[i+1]-ageStratification[i]));
                eta.push_back( 1.0/(365.0*(ageStratification[i+1] - ageStratification[i])) );
                modelIncidencePerTime.push_back(0);
            }
            modelIncidencePerTime.push_back(0);
            populationPerAgeGroup.push_back(totPopulation - (dailyBirthRate*365)*ageStratification[A-1]);
            eta.push_back(dailyBirthRate/(totPopulation - (dailyBirthRate*365)*ageStratification[A-1]) );

            dt = 1;
            currentODETime = 0;
            dayNoAfterBurn = 0;
            valueLogLikelihood = 0;  
    }
    
    vector<double > pVHR;
    vector<double > pHR;
    vector<double > pLR;
    
    // Related to parameter values
    NumericVector parameterValues;
    
    void ParameterValuesforODE(VectorXd currentParamValues){

        this->ep_t.clear();
        this->pA.clear();
        
        NumericVector parameterValuesTemp(16);
        
      /*   parameterValuesTemp["xi"] = currentParamValues(0);
         parameterValuesTemp["si"] = currentParamValues(1);
         parameterValuesTemp["ga0"] = currentParamValues(2);
         parameterValuesTemp["g1"] = currentParamValues(3);
         parameterValuesTemp["g2"] = currentParamValues(4);
         parameterValuesTemp["om"] = currentParamValues(5);
         parameterValuesTemp["pA1"] = currentParamValues(6);
         parameterValuesTemp["pA2"] = currentParamValues(7);
         parameterValuesTemp["pA3"] = currentParamValues(8);
         parameterValuesTemp["pA4"] = currentParamValues(9);*/
         parameterValuesTemp["alpha_i"] = currentParamValues(0);
         /*
         parameterValuesTemp["d1"] = currentParamValues(11);
         parameterValuesTemp["d2"] = currentParamValues(12);
         parameterValuesTemp["d3"] = currentParamValues(13);*/
         parameterValuesTemp["phi"] = currentParamValues(1);
         parameterValuesTemp["qp"] = currentParamValues(2);
         parameterValuesTemp["qc"] = currentParamValues(3);
         parameterValuesTemp["b1"] = currentParamValues(4);
         parameterValuesTemp["psi"] = currentParamValues(5);
         parameterValuesTemp["ep1"] = currentParamValues(6);
         parameterValuesTemp["ep2"] = currentParamValues(7);
         parameterValuesTemp["ep3"] = currentParamValues(8);
         parameterValuesTemp["ep4"] = currentParamValues(9);
         parameterValuesTemp["ep5"] = currentParamValues(10);
         parameterValuesTemp["ep6"] = currentParamValues(11);
         parameterValuesTemp["ep7"] = currentParamValues(12);
         parameterValuesTemp["ep8"] = currentParamValues(13);
         parameterValuesTemp["ep9"] = currentParamValues(14);
         parameterValuesTemp["I1"] = currentParamValues(15);
         parameterValuesTemp["I2"] = currentParamValues(16);
        this->parameterValues = parameterValuesTemp;
        
        // Define ep_t;
                
        for (int a = 0; a < 4; a++)
            this->ep_t.push_back(parameterValues["ep1"]);
        for (int a = 4; a < 7; a++)
            this->ep_t.push_back(parameterValues["ep2"]);
        for (int a = 7; a < 12; a++)
            this->ep_t.push_back(parameterValues["ep3"]);
        for (int a = 12; a < 13; a++)
            this->ep_t.push_back(parameterValues["ep4"]);
        for (int a = 13; a < 16; a++)
            this->ep_t.push_back(parameterValues["ep5"]);
        for (int a = 16; a < 22; a++)
            this->ep_t.push_back(parameterValues["ep6"]);
        for (int a = 22; a < 23; a++)
            this->ep_t.push_back(parameterValues["ep7"]);
        for (int a = 23; a < 24; a++)
            this->ep_t.push_back(parameterValues["ep8"]);
        for (int a = 24; a < 25; a++)
            this->ep_t.push_back(parameterValues["ep9"]);

        /*
        for (int a = 0; a < 16; a++)
            this->ep_t.push_back(exp(parameterValues["c5ep1"] + a*parameterValues["c5ep2"]));
        for (int a = 16; a < 23; a++)
            this->ep_t.push_back(parameterValues["ep5"]);
        for (int a = 23; a < 25; a++)
            this->ep_t.push_back(parameterValues["ep6"]);*/

        for (int a = 0; a < 12; a++)
            this->pA.push_back(0.0916);
        for (int a = 12; a < 16; a++)
            this->pA.push_back(0.163);
        for (int a = 16; a < 18; a++)
            this->pA.push_back(0.516);
        for (int a = 18; a < 25; a++)
            this->pA.push_back(0.753);
        // define pA;
        /*
        for (int a = 0; a < 12; a++)
            this->pA.push_back(parameterValues["pA1"]);
        for (int a = 12; a < 16; a++)
            this->pA.push_back(parameterValues["pA2"]);
        for (int a = 16; a < 18; a++)
            this->pA.push_back(parameterValues["pA3"]);
        for (int a = 18; a < 25; a++)
            this->pA.push_back(parameterValues["pA4"]);*/
    }
    
    vector<double > initial_M()
    {
       // double xi = 1.0/this->parameterValues["xi"];
        double xi = 1.0/60.0;
        boost::math::exponential_distribution <> exp ( xi );
        
        vector<double > init_con;
        for (int i = 0; i< this->A - 1; i++)
        {
            double init_con_temp = (cdf(exp, 365*ageStratification[i+1])-cdf(exp, 365*ageStratification[i]))/((365*ageStratification[i+1]-365*ageStratification[i])*xi);
            init_con.push_back(init_con_temp*populationPerAgeGroup[i]);
        }
        init_con.push_back( (cdf(exp, 365*90)-cdf(exp, 365*ageStratification[A-1]))/((365*90-365*ageStratification[A-1])*xi)*populationPerAgeGroup[A-1]);
        return init_con;
    }

    double poisson_cdf(double l, double a, double x){
      if( l == 0.0 || a == 0.0){
        boost::math::poisson_distribution<> p(0.000001); return cdf(p,x);
      }
      else{
        boost::math::poisson_distribution<> p(l*a); return cdf(p,x);
      }
    }

    vector< double >  initialProportionExposure(double l, double a1, double a2){
        vector< double >  prop(this->A);
        prop[0] = abs(poisson_cdf(l,a2,0)-poisson_cdf(l,a1,0))/((a2-a1)*l);
        prop[1] = abs(poisson_cdf(l,a2,1)-poisson_cdf(l,a1,1))/((a2-a1)*l);
        prop[2] = abs(poisson_cdf(l,a2,2)-poisson_cdf(l,a1,2))/((a2-a1)*l);
        prop[3] = 1 - (prop[2]+prop[1]+prop[0]);
        return prop;
    }

    vector<double > generateInitialStates(double cov_c)
    {
        vector<double > initialStates;
        vector<double > populationMatPro = initial_M();
        
        double s_prop, r_prop;
        double pI1; double pI2; double pI3; double pI4;
        
        double I1 = this->parameterValues["I1"];
        double I2 = this->parameterValues["I2"];
        double I3 = 0.5;
        
    /*   double si = 1.0/parameterValues["si"];
        double g0 = 1.0/parameterValues["ga0"];
        double g1 = 1.0/((parameterValues["ga0"])*(parameterValues["g1"]));
        double g2 = 1.0/((parameterValues["ga0"])*(parameterValues["g1"])*(parameterValues["g2"]));
        double d1 = parameterValues["d1"];
        double d2 = parameterValues["d1"]*parameterValues["d2"];
        double d3 = parameterValues["d1"]*parameterValues["d2"]*parameterValues["d3"];*/

        double si = 1.0/4.98;
        double g0 = 1.0/6.16;
        double g1 = 1.0/(6.16 * 0.87);
        double g2 = 1.0/(6.16 * 0.87 * 0.79);
        double g3 = g2;

        double d1 = 0.89;
        double d2 = 0.89 * 0.81;
        double d3 = 0.89 * 0.81 * 0.33;

        double a1, a2;
        
        for(int a = 0; a < A; a++)
        {
            if (a < A-1){
                a1 = this->ageStratification(a); a2 =  this->ageStratification(a+1);
            }
            else{
                a1 =  this->ageStratification(a); a2 = 90;
            }
            
            vector<double > propEachExposureGroup = initialProportionExposure(I3, a1, a2);
            pI1 = propEachExposureGroup[0]; pI2 = propEachExposureGroup[1]; pI3 = propEachExposureGroup[2]; pI4 = propEachExposureGroup[3];
            /** Get ratio of number of number of each person in each compoartment group **/
            double populationNoMatPro = populationPerAgeGroup[a] - populationMatPro[a];
            
            for (int s = 0; s < 6; s++)
            {
                if (a < 12)
                {
                    if (s == 0 || s == 3)
                        s_prop = (0.0);
                    else if (s == 1 || s == 4)
                        s_prop = (cov_c);
                    else if (s == 2 || s == 5)
                        s_prop = (1-cov_c);
                    else
                        cout << "OOPS" << '\n';
                }
                else
                {
                    if (s == 0 || s == 3)
                        s_prop = p_mat[a]*(1-cov_c);
                    else if (s == 1 || s == 4)
                        s_prop = p_mat[a]*(cov_c);
                    else if (s == 2 || s == 5)
                        s_prop = 1-p_mat[a];
                    
                    else
                        cout << "OOPS" << '\n';
                }
                for (int r = 0; r < 3; r++)
                {
                    if (r == 0)
                        r_prop = pVHR[a];
                    else if (r == 1)
                        r_prop = pHR[a];
                    else if (r == 2)
                        r_prop = pLR[a];
                    else
                        cout << "ERROR" << endl;
                    
                    initialStates.push_back(r_prop*s_prop*populationMatPro[a]); // Number  in M group
                    
                    initialStates.push_back(r_prop*s_prop*pI1*populationNoMatPro*(1.0 - I1)*(1.0-I2));      //Sus
                    initialStates.push_back(r_prop*s_prop*pI1*populationNoMatPro*I1*si/(si+g0));       //Exp
                    initialStates.push_back(r_prop*s_prop*pI1*populationNoMatPro*I1*g0/(si+g0)*pA[a]);
                    initialStates.push_back(r_prop*s_prop*pI1*populationNoMatPro*I1*g0/(si+g0)*(1-pA[a]));   //Inf S
                    initialStates.push_back(r_prop*s_prop*pI1*populationNoMatPro*(1.0 - I1)*I2);      //Rec
                    
                    initialStates.push_back(r_prop*s_prop*pI2*populationNoMatPro*(1.0 - d1*I1)*(1.0-I2));      //Sus
                    initialStates.push_back(r_prop*s_prop*pI2*populationNoMatPro*d1*I1*si/(si+g1));       //Exp
                    initialStates.push_back(r_prop*s_prop*pI2*populationNoMatPro*d1*I1*g1/(si+g1)*pA[a]);
                    initialStates.push_back(r_prop*s_prop*pI2*populationNoMatPro*d1*I1*g1/(si+g1)*(1-pA[a]));   //Inf S
                    initialStates.push_back(r_prop*s_prop*pI2*populationNoMatPro*(1.0 - d1*I1)*I2);      //Rec
                    
                    initialStates.push_back(r_prop*s_prop*pI3*populationNoMatPro*(1.0 - d2*I1)*(1.0-I2));      //S/Sus
                    initialStates.push_back(r_prop*s_prop*pI3*populationNoMatPro*d2*I1*si/(si+g2));       //Exp
                    initialStates.push_back(r_prop*s_prop*pI3*populationNoMatPro*d2*I1*g2/(si+g2)*pA[a]);
                    initialStates.push_back(r_prop*s_prop*pI3*populationNoMatPro*d2*I1*g2/(si+g2)*(1-pA[a]));   //Inf S
                    initialStates.push_back(r_prop*s_prop*pI3*populationNoMatPro*(1.0 - d2*I1)*I2);      //ec
                    
                    initialStates.push_back(r_prop*s_prop*pI4*populationNoMatPro*(1.0 - d3*I1)*(1.0-I2));        //Sus
                    initialStates.push_back(r_prop*s_prop*pI4*populationNoMatPro*d3*I1*si/(si+g2));       //Exp
                    initialStates.push_back(r_prop*s_prop*pI4*populationNoMatPro*d3*I1*g2/(si+g2)*pA[a]);
                    initialStates.push_back(r_prop*s_prop*pI4*populationNoMatPro*d3*I1*g2/(si+g2)*(1-pA[a]));   //Inf S
                    initialStates.push_back(r_prop*s_prop*pI4*populationNoMatPro*(1.0 - d3*I1)*I2);     //Rec
                    
                    initialStates.push_back(0);     //Rec
                    initialStates.push_back(0);     //Rec
                    initialStates.push_back(0);     //Rec

                }
            }
            for (int j = 0; j < 23; j++)
                initialStates.push_back(0.0);      //inf exposure 1
        }
        return initialStates;
    }
    
    void getWeeklyIncidence(vector<double> &x0, NumericMatrix &sampleWeeklyIncidence, NumericMatrix &no_doses, bool epFlag)
    {
      if (this->dayNoAfterBurn == 0){
        for (int a = 0; a < this->A; a++) {
            x0[455*a + 72*6 + 1] = 0;
            x0[455*a + 72*6 + 2] = 0;
            x0[455*a + 72*6 + 3] = 0;
            x0[455*a + 72*6 + 4] = 0;
            for (int j = 0; j < 9; j++) {
                x0[455*a + 72*6 + 8 + j] = 0.0; //Incidence at t_d = 0;
            }
        }
      }
        if (this->dayNoAfterBurn%7 == 0 && this->dayNoAfterBurn > 0)
        {
            for (int a = 0; a < this->A; a++) {
                for (int j = 0; j < 9; j++)
                {
                    if (epFlag)
                      sampleWeeklyIncidence(this->weekNo, 9*a + j) = x0[455*a + 72*6 + 8 + j]*this->ep_t[a];
                    else
                      sampleWeeklyIncidence(this->weekNo, 9*a + j) = x0[455*a + 72*6 + 8 + j]; //Incidence at t_d = 7;
                  
                    x0[455*a + 72*6 + 8 + j] = 0.0;
                }
            }
            for (int a = 0; a < this->A; a++) {
                no_doses(this->weekNo, 0) += (x0[455*a + 72*6 + 1]); // pal
                no_doses(this->weekNo, 1) += (x0[455*a + 72*6 + 2]); // mab
                no_doses(this->weekNo, 2) += (x0[455*a + 72*6 + 3]); // lav
                no_doses(this->weekNo, 3) += (x0[455*a + 72*6 + 4]); // mat
                
                x0[455*a + 72*6 + 1] = 0;
                x0[455*a + 72*6 + 2] = 0;
                x0[455*a + 72*6 + 3] = 0;
                x0[455*a + 72*6 + 4] = 0;
            }

            this->weekNo++;
        }
        this->dayNoAfterBurn++;
    }

    void getMonthlyIncidence(vector<double> &x0, NumericMatrix &sampleMonthlyIncidence, NumericMatrix &no_doses, bool epFlag)
    {
      if (this->dayNoAfterBurn == 0){
        for (int a = 0; a < this->A; a++) {
            x0[455*a + 72*6 + 1] = 0;
            x0[455*a + 72*6 + 2] = 0;
            x0[455*a + 72*6 + 3] = 0;
            x0[455*a + 72*6 + 4] = 0;
            for (int j = 0; j < 9; j++) {
                x0[455*a + 72*6 + 8 + j] = 0.0; //Incidence at t_d = 0;
            }
        }
      }
        if (this->dayNoAfterBurn%30 == 0 && this->dayNoAfterBurn > 0)
        {
            for (int a = 0; a < this->A; a++) {
                for (int j = 0; j < 9; j++)
                {
                    if (epFlag)
                      sampleMonthlyIncidence(this->monthNo, 9*a + j) = x0[455*a + 72*6 + 8 + j]*this->ep_t[a];
                    else
                      sampleMonthlyIncidence(this->monthNo, 9*a + j) = x0[455*a + 72*6 + 8 + j]; //Incidence at t_d = 7;
                  
                    x0[455*a + 72*6 + 8 + j] = 0.0;
                }
            }
            for (int a = 0; a < this->A; a++) {
                no_doses(this->monthNo, 0 + 4 * (a)) = (x0[455*a + 72*6 + 1]); // pal
                no_doses(this->monthNo, 1 + 4 * (a)) = (x0[455*a + 72*6 + 2]); // mab
                no_doses(this->monthNo, 2 + 4 * (a)) = (x0[455*a + 72*6 + 3]); // lav
                no_doses(this->monthNo, 3 + 4 * (a)) = (x0[455*a + 72*6 + 4]); // mat
                
                x0[455*a + 72*6 + 1] = 0;
                x0[455*a + 72*6 + 2] = 0;
                x0[455*a + 72*6 + 3] = 0;
                x0[455*a + 72*6 + 4] = 0;
            }
            this->monthNo++;
        }
        this->dayNoAfterBurn++;
    }
};


class ODE_desc
{
    RunInterventions* finELL;
  
public:
    double M, S0, S1, S2, S3, SN, E0, E1, E2, E3, A0, A1, A2, A3, I0, I1, I2, I3, R0, R1, R2, R3, M_, S0_, S1_, S2_, S3_, E0_, E1_, E2_, E3_, A0_, A1_, A2_, A3_, I0_, I1_, I2_, I3_, R0_, R1_, R2_, R3_, N;
    double Mv, S0v, S1v, S2v, S3v, E0v, E1v, E2v, E3v, A0v, A1v, A2v, A3v, I0v, I1v, I2v, I3v, R0v, R1v, R2v, R3v, Mv_, S0v_, S1v_, S2v_, S3v_, E0v_, E1v_, E2v_, E3v_, A0v_, A1v_, A2v_, A3v_, I0v_, I1v_, I2v_, I3v_, R0v_, R1v_, R2v_, R3v_;
    double xi, si, ga0, ga1, ga2, ga3, d1, d2, d3, a1, a2, a3, alpha_i, rho, om, b1, qp, qc, psi, phi, beta;
    int A;
    double p_vul;
    VectorXd I_n;     VectorXd I_c;     VectorXd I_p;
    VectorXd I_n_v;   VectorXd I_c_v;   VectorXd I_p_v;
    
    VectorXd N_tot;
    VectorXd N_tot_n;         VectorXd N_tot_c;         VectorXd N_tot_p;
    VectorXd N_tot_n_v;       VectorXd N_tot_c_v;       VectorXd N_tot_p_v;
    double N_tot_n_v_inv, N_tot_p_v_inv, N_tot_c_v_inv, N_tot_n_inv, N_tot_p_inv, N_tot_c_inv;
    
    VectorXd prop_n;  VectorXd prop_c;  VectorXd prop_p;
    VectorXd prop_nv; VectorXd prop_cv; VectorXd prop_pv;
    VectorXd prop_empty;
    
    VectorXd PS;
    
    double I_temp_p, I_temp_c, I_temp_n, I_temp_p_v, I_temp_c_v, I_temp_n_v;
    double x_tot, x_tot_1, x_tot_2;
    double xi_b;
    
    vector<double > pVHR;
    vector<double > pHR;
    vector<double > pLR;
    vector<double > p_mat;

    double dailyBirthRate;
    NumericMatrix contactMatrixPhy, contactMatrixCon;
    vector< double >  populationPerAgeGroup, eta, pA;
    vector<double > VHR_g = {0.653012, 0.291842, 1., 1., 1., 0.0321774, 1., 1., 0,0,0,0,0,0,0,0,0,0,0};
    vector<double > LR_g = {1.00126, 1.00168, 1., 1., 1., 1.00067, 1., 1., 1.00002, 1, 1, 1, 1, 1, 1, 1, 1};
    
    MatrixXd cal_pal, cal_mAB_VHR, cal_mAB_HR, cal_mAB_LR, cal_LAV_HR, cal_LAV_LR, cal_mat_LR;
    MatrixXd cal_pal_t, cal_mAB_VHR_t, cal_mAB_HR_t, cal_mAB_LR_t, cal_LAV_HR_t, cal_LAV_LR_t, cal_mat_LR_t;
    MatrixXd cal_pal_dose, cal_mAB_VHR_dose, cal_mAB_HR_dose, cal_mAB_LR_dose, cal_LAV_HR_dose, cal_LAV_LR_dose, vac_cal_dose;
    List vac_calendar, vac_info, vac_dose;
    double cov_c, om_mab, xi_boost;
    bool direct;

    double u18p;
    double u19p;
    double u20p;
    
  ODE_desc(RunInterventions* finELL_t, List vac_calendar_t, List vac_dose_t, List vac_info_t, double cov_c_t): finELL(finELL_t), vac_calendar(vac_calendar_t), vac_dose(vac_dose_t), vac_info(vac_info_t), cov_c(cov_c_t){

    //  xi = 1.0/finELL->parameterValues["xi"];
    //  si = 1.0/finELL->parameterValues["si"];
    //  ga0 = 1.0/(finELL->parameterValues["ga0"]);
    //   ga1 = 1.0/(finELL->parameterValues["ga0"]*finELL->parameterValues["g1"]);
    //  ga2 = 1.0/(finELL->parameterValues["ga0"]*finELL->parameterValues["g1"]*finELL->parameterValues["g2"]);
    //  ga3 = ga2;
    //  om = 1.0/finELL->parameterValues["om"];
        xi = 1.0/60.0;
        si = 1.0/4.98;
        ga0 = 1.0/6.16;
        ga1 = 1.0/(6.16 * 0.87);
        ga2 = 1.0/(6.16 * 0.87 * 0.79);
        ga3 = ga2;
        om = 1.0/358.4;

        rho = 1.0;
        //alpha_i = 0.634;
        alpha_i = finELL->parameterValues["alpha_i"];
    // d1 = finELL->parameterValues["d1"];
    // d2 = finELL->parameterValues["d1"]*finELL->parameterValues["d2"];
    //  d3 = finELL->parameterValues["d1"]*finELL->parameterValues["d2"]*finELL->parameterValues["d3"];
        d1 = 0.89;
        d2 = 0.89 * 0.81;
        d3 = 0.89 * 0.81 * 0.33;
        a1 = 1.0, a2 = 1.0, a3 = 1.0;

        phi = finELL->parameterValues["phi"];
        qp = finELL->parameterValues["qp"];
        qc = finELL->parameterValues["qc"];
        b1 = finELL->parameterValues["b1"];
        psi = finELL->parameterValues["psi"];

        A = finELL->A;
        eta = finELL->eta;
      
        u18p = finELL->u_p[0];
        u19p = finELL->u_p[1];
        u20p = finELL->u_p[2];

        populationPerAgeGroup = finELL->populationPerAgeGroup;
        pVHR = finELL->pVHR;
        pHR = finELL->pHR;
        pLR = finELL->pLR;
        p_mat = finELL->p_mat;

        dailyBirthRate = finELL->dailyBirthRate;
        pA = finELL->pA;
      
        I_n = VectorXd::Zero(A);
        I_c = VectorXd::Zero(A);
        I_p = VectorXd::Zero(A);
        I_n_v = VectorXd::Zero(A);
        I_c_v = VectorXd::Zero(A);
        I_p_v = VectorXd::Zero(A);
      
        N_tot = VectorXd::Zero(A);
        N_tot_n = VectorXd::Zero(A);
        N_tot_c = VectorXd::Zero(A);
        N_tot_p = VectorXd::Zero(A);
        N_tot_n_v = VectorXd::Zero(A);
        N_tot_c_v = VectorXd::Zero(A);
        N_tot_p_v = VectorXd::Zero(A);
      
        prop_n = VectorXd::Zero(A);   prop_c = VectorXd::Zero(A);   prop_p = VectorXd::Zero(A);
        prop_nv = VectorXd::Zero(A);  prop_cv = VectorXd::Zero(A);  prop_pv = VectorXd::Zero(A);
        prop_empty = VectorXd::Zero(A);
      
        PS = VectorXd::Zero(A - 1);

        // This needs sorting out somehow
        cal_pal_t = vac_calendar["pal"];
        cal_mAB_VHR_t = vac_calendar["mAB_VHR"];
        cal_mAB_HR_t = vac_calendar["mAB_HR"];
        cal_mAB_LR_t = vac_calendar["mAB_LR"];
        cal_LAV_HR_t = vac_calendar["LAV_HR"];
        cal_LAV_LR_t = vac_calendar["LAV_LR"];
        cal_mat_LR_t = vac_calendar["mat_LR"];
      
        cal_pal_dose = vac_dose["pal"];
        cal_mAB_VHR_dose = vac_dose["mAB_VHR"];
        cal_mAB_HR_dose = vac_dose["mAB_HR"];
        cal_mAB_LR_dose = vac_dose["mAB_LR"];
        cal_LAV_HR_dose = vac_dose["LAV_HR"];
        cal_LAV_LR_dose = vac_dose["LAV_LR"];
        vac_cal_dose = vac_dose["mat_LR"];
      
        //cov_c = cov_c_t; 
        direct = vac_info["direct"];
        om_mab = vac_info["om_mab"];
        xi_boost = vac_info["xi_boost"];
   };
  
    
  
  void operator() (  vector< double >  &x , vector< double >  &dxdt , const double  t )
  {
      int ag = 455;
      int sg = 72;
      int rg = 24;
      /*/////////////*/
      /* CALENDARS /*/
      /*/////////////*/
    //  Rcpp::Rcout << "Defining the Calendars" << std::endl;
      MatrixXd vac_cal_vhr; MatrixXd vac_cal;
      if (t < 365)
      {
          cal_pal = MatrixXd::Zero(365, A);
          vac_cal = MatrixXd::Zero(365, A);
          cal_mAB_VHR = MatrixXd::Zero(365, A);
          cal_mAB_HR = MatrixXd::Zero(365, A);
          cal_mAB_LR = MatrixXd::Zero(365, A);
          cal_LAV_HR = MatrixXd::Zero(365, A);
          cal_LAV_LR = MatrixXd::Zero(365, A);
      }
      else{
          cal_pal = cal_pal_t;
          vac_cal = cal_mat_LR_t;
          cal_mAB_VHR = cal_mAB_VHR_t;
          cal_mAB_HR = cal_mAB_HR_t;
          cal_mAB_LR = cal_mAB_LR_t;
          cal_LAV_HR = cal_LAV_HR_t;
          cal_LAV_LR = cal_LAV_LR_t;
      }
      
      /*/////////////*/
      /* POPULATIONS /*/
      /*/////////////*/
      //Find sum of every group
   //   Rcpp::Rcout << "Info on popoulations" << std::endl;

      double tot_coc = 0, tot_coc_v = 0, tot_inf = 0, tot_inf_v = 0;
      
      for(int a = 0; a < A; a++)
      {
          double r_prop = 0;
          double tot_temp_n = 0, tot_temp_c = 0, tot_temp_p = 0, tot_temp_nv = 0, tot_temp_cv = 0, tot_temp_pv = 0;
          double tot_temp_vhr = 0, tot_temp_vhrv = 0, tot_temp_hr = 0, tot_temp_hrv = 0, tot_temp_lr = 0, tot_temp_lrv = 0;

          for (int j = 0; j < 24; j++){
              for (int r = 0; r < 3; r++){
                  tot_temp_p +=  x[a*ag+0*sg+r*rg+j];        tot_temp_c += x[a*ag+1*sg+r*rg+j];            tot_temp_n += x[a*ag+2*sg+r*rg+j];
                  tot_temp_pv += x[a*ag+3*sg+r*rg+j];       tot_temp_cv += x[a*ag+4*sg+r*rg+j];           tot_temp_nv += x[a*ag+5*sg+r*rg+j];
              }
          }
          
          for (int j = 0; j < 24; j++){
              for (int s = 0; s < 3; s++){
                  tot_temp_vhr +=  x[a*ag+s*sg+0*rg+j];        tot_temp_hr += x[a*ag+s*sg+1*rg+j];            tot_temp_lr += x[a*ag+s*sg+2*rg+j];
                  tot_temp_vhrv += x[a*ag+(s+3)*sg+0*rg+j];       tot_temp_hrv += x[a*ag+(s+3)*sg+1*rg+j];           tot_temp_lrv += x[a*ag+(s+3)*sg+2*rg+j];
              }
          }
          
          if (a < 12){
              tot_coc += tot_temp_c;                  tot_coc_v += tot_temp_cv;
              tot_inf += tot_temp_c + tot_temp_n;     tot_inf_v += tot_temp_cv + tot_temp_nv;
          }
      
          N_tot_n[a] = tot_temp_n;    N_tot_c[a] = tot_temp_c;    N_tot_p[a] = tot_temp_p;
          N_tot_n_v[a] = tot_temp_nv; N_tot_c_v[a] = tot_temp_cv; N_tot_p_v[a] = tot_temp_pv;

      }
      
      /*/////////////*/
      /* 2. Force of Infection /*/
      /*/////////////*/
    //  Rcpp::Rcout << "Info on FOI" << std::endl;

      double phi_c = cov_c;
      double phi_c_v = cov_c;
      
      //  2.1 Get contact matrices
      vector<vector< double > > cnt_matrix_cwc_p = finELL->get_cwc(phi_c,'p');
      vector<vector< double > > cnt_matrix_cwp_p = finELL->get_cwp(phi_c,'p');
      vector<vector< double > > cnt_matrix_cwn_p = finELL->get_cwn(phi_c,'p');
      vector<vector< double > > cnt_matrix_pwc_p = finELL->get_pwc(phi_c,'p');
      vector<vector< double > > cnt_matrix_pwp_p = finELL->get_pwp(phi_c,'p');
      vector<vector< double > > cnt_matrix_pwn_p = finELL->get_pwn(phi_c,'p');
      vector<vector< double > > cnt_matrix_nwc_p = finELL->get_nwc(phi_c,'p');
      vector<vector< double > > cnt_matrix_nwp_p = finELL->get_nwp(phi_c,'p');
      vector<vector< double > > cnt_matrix_nwn_p = finELL->get_nwn(phi_c,'p');
      
      vector<vector< double > > cnt_matrix_cwc_c = finELL->get_cwc(phi_c,'c');
      vector<vector< double > > cnt_matrix_cwp_c = finELL->get_cwp(phi_c,'c');
      vector<vector< double > > cnt_matrix_cwn_c = finELL->get_cwn(phi_c,'c');
      vector<vector< double > > cnt_matrix_pwc_c = finELL->get_pwc(phi_c,'c');
      vector<vector< double > > cnt_matrix_pwp_c = finELL->get_pwp(phi_c,'c');
      vector<vector< double > > cnt_matrix_pwn_c = finELL->get_pwn(phi_c,'c');
      vector<vector< double > > cnt_matrix_nwc_c = finELL->get_nwc(phi_c,'c');
      vector<vector< double > > cnt_matrix_nwp_c = finELL->get_nwp(phi_c,'c');
      vector<vector< double > > cnt_matrix_nwn_c = finELL->get_nwn(phi_c,'c');

      vector<vector< double > > cnt_matrix_cwc_p_v = finELL->get_cwc(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_cwp_p_v = finELL->get_cwp(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_cwn_p_v = finELL->get_cwn(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_pwc_p_v = finELL->get_pwc(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_pwp_p_v = finELL->get_pwp(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_pwn_p_v = finELL->get_pwn(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_nwc_p_v = finELL->get_nwc(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_nwp_p_v = finELL->get_nwp(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_nwn_p_v = finELL->get_nwn(phi_c_v,'p');
      
      vector<vector< double > > cnt_matrix_cwc_c_v = finELL->get_cwc(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_cwp_c_v = finELL->get_cwp(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_cwn_c_v = finELL->get_cwn(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_pwc_c_v = finELL->get_pwc(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_pwp_c_v = finELL->get_pwp(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_pwn_c_v = finELL->get_pwn(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_nwc_c_v = finELL->get_nwc(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_nwp_c_v = finELL->get_nwp(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_nwn_c_v = finELL->get_nwn(phi_c_v,'c');
      
      //  2.2 Seasonal forcing
      //beta = (1 + b1*cos((t/365.0-phi)*2*PI));
      int t1 = (int)t%365;
      beta = 1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi)));

      /*/////////////*/
      /*  3. DYNAMIC MATERNAL PROTECTION /*/
      /*/////////////*/
  //    Rcpp::Rcout << "Dynamic maternal protection" << std::endl;

     // num_vec pos_wcb_ag = {18,21};
      double sum_wcb = 0;
      double sum_wcb_v = 0;
      double CB2_temp, CB2_temp_v;
      CB2_temp = CB2_temp_v = 0;
      for(int a = 18; a < 21 ; a++)
      {
          for (int r = 0; r < 3; r++)
          {
              CB2_temp += (x[a*ag+2*sg+r*rg+1] + x[a*ag+2*sg+r*rg+6] + x[a*ag+2*sg+r*rg+11] + x[a*ag+2*sg+r*rg+16] + x[a*ag+2*sg+r*rg+2] + x[a*ag+2*sg+r*rg+7] + x[a*ag+2*sg+r*rg+12] + x[a*ag+2*sg+r*rg+17]);
              CB2_temp_v += (x[a*ag+5*sg+r*rg+1] + x[a*ag+5*sg+r*rg+6] + x[a*ag+5*sg+r*rg+11] + x[a*ag+5*sg+r*rg+16] + x[a*ag+5*sg+r*rg+2] + x[a*ag+5*sg+r*rg+7] + x[a*ag+5*sg+r*rg+12] + x[a*ag+5*sg+r*rg+17]);
          }
      }
      sum_wcb = CB2_temp/((double)N_tot_n[18]+(double)N_tot_n[19]+(double)N_tot_n[20]);
      sum_wcb_v = (CB2_temp_v)/((double)N_tot_n[18]+(double)N_tot_n[19]+(double)N_tot_n[20]);
      
      double tot_C = N_tot_c[18] + N_tot_c[19] + N_tot_c[20];
      double tot_C_inv;
      if (tot_C < 1.0)
          tot_C_inv = 0;
      else
          tot_C_inv = 1.0/tot_C;
      /*/////////////*/
      /* 4. ODES /*/
      /*/////////////*/
      //4.1 Age groups
      double protectpal;
      double protectmabs;
      double protectLAV;
      double protectmat;

    //  Rcpp::Rcout << "ODEs" << std::endl;

      for(int a = 0; a < A; a++)
      {
          double r_prop = 0;
          protectpal = protectmabs = protectLAV = protectmat = 0;
          I_temp_c = 0.0;    I_temp_p = 0.0; I_temp_n = 0.0;
      //    Rcpp::Rcout << "Defining the FOI for age group: " << a << std::endl;

          for (int k = 0; k < A; k++)
          {
              
              if (N_tot_n[k] < 0.1){N_tot_n_inv = 0;}else{N_tot_n_inv=1.0/N_tot_n[k];}
              if (N_tot_c[k] < 0.1){N_tot_c_inv = 0;}else{N_tot_c_inv=1.0/N_tot_c[k];}
              if (N_tot_p[k] < 0.1){N_tot_p_inv = 0;}else{N_tot_p_inv=1.0/N_tot_p[k];}
              
              for (int r = 0; r < 3; r++)
              {
                  I_temp_p += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_pwp_p[a][k]+qc*cnt_matrix_pwp_c[a][k]))*N_tot_p_inv +(x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_pwc_p[a][k]+qc*cnt_matrix_pwc_c[a][k]))*N_tot_c_inv +(x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_pwn_p[a][k]+qc*cnt_matrix_pwn_c[a][k]))*N_tot_n_inv;
                  
                  I_temp_c += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_cwp_p[a][k]+qc*cnt_matrix_cwp_c[a][k]))*N_tot_p_inv +(x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_cwc_p[a][k]+qc*cnt_matrix_cwc_c[a][k]))*N_tot_c_inv +(x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_cwn_p[a][k]+qc*cnt_matrix_cwn_c[a][k]))*N_tot_n_inv;
                  
                  I_temp_n += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_nwp_p[a][k]+qc*cnt_matrix_nwp_c[a][k]))*N_tot_p_inv +(x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_nwc_p[a][k]+qc*cnt_matrix_nwc_c[a][k]))*N_tot_c_inv +(x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_nwn_p[a][k]+qc*cnt_matrix_nwn_c[a][k]))*N_tot_n_inv;
              }

          }

          
          I_temp_c_v = 0.0;  I_temp_p_v = 0.0; I_temp_n_v = 0.0;
          for (int k = 0; k < A ; k++)
          {
              if (N_tot_n_v[k] < 1){N_tot_n_v_inv = 0;}else{N_tot_n_v_inv=1.0/N_tot_n_v[k];}
              if (N_tot_c_v[k] < 1){N_tot_c_v_inv = 0;}else{N_tot_c_v_inv=1.0/N_tot_c_v[k];}
              if (N_tot_p_v[k] < 1){N_tot_p_v_inv = 0;}else{N_tot_p_v_inv=1.0/N_tot_p_v[k];}
              
              for (int r = 0; r < 3; r++)
              {
                  I_temp_p_v += (x[k*ag+3*sg+r*rg+3]*alpha_i+x[k*ag+3*sg+r*rg+4]+a1*(x[k*ag+3*sg+r*rg+8]*alpha_i+x[k*ag+3*sg+r*rg+9])+a2*(x[k*ag+3*sg+r*rg+13]*alpha_i+x[k*ag+3*sg+r*rg+14])+a3*(x[k*ag+3*sg+r*rg+18]*alpha_i+x[k*ag+3*sg+r*rg+19]))*(qp*(cnt_matrix_pwp_p_v[a][k]+qc*cnt_matrix_pwp_c_v[a][k]))*N_tot_p_v_inv + (x[k*ag+4*sg+r*rg+3]*alpha_i+x[k*ag+4*sg+r*rg+4]+a1*(x[k*ag+4*sg+r*rg+8]*alpha_i+x[k*ag+4*sg+r*rg+9])+a2*(x[k*ag+4*sg+r*rg+13]*alpha_i+x[k*ag+4*sg+r*rg+14])+a3*(x[k*ag+4*sg+r*rg+18]*alpha_i+x[k*ag+4*sg+r*rg+19]))*(qp*(cnt_matrix_pwc_p_v[a][k]+qc*cnt_matrix_pwc_c_v[a][k]))*N_tot_c_v_inv  + (x[k*ag+5*sg+r*rg+3]*alpha_i+x[k*ag+5*sg+r*rg+4]+a1*(x[k*ag+5*sg+r*rg+8]*alpha_i+x[k*ag+5*sg+r*rg+9])+a2*(x[k*ag+5*sg+r*rg+13]*alpha_i+x[k*ag+5*sg+r*rg+14])+a3*(x[k*ag+5*sg+r*rg+18]*alpha_i+x[k*ag+5*sg+r*rg+19]))*(qp*(cnt_matrix_pwn_p_v[a][k]+qc*cnt_matrix_pwn_c_v[a][k]))*N_tot_n_v_inv;
                  
                  I_temp_c_v += (x[k*ag+3*sg+r*rg+3]*alpha_i+x[k*ag+3*sg+r*rg+4]+a1*(x[k*ag+3*sg+r*rg+8]*alpha_i+x[k*ag+3*sg+r*rg+9])+a2*(x[k*ag+3*sg+r*rg+13]*alpha_i+x[k*ag+3*sg+r*rg+14])+a3*(x[k*ag+3*sg+r*rg+18]*alpha_i+x[k*ag+3*sg+r*rg+19]))*(qp*(cnt_matrix_cwp_p_v[a][k]+qc*cnt_matrix_cwp_c_v[a][k]))*N_tot_p_v_inv +(x[k*ag+4*sg+r*rg+3]*alpha_i+x[k*ag+4*sg+r*rg+4]+a1*(x[k*ag+4*sg+r*rg+8]*alpha_i+x[k*ag+4*sg+r*rg+9])+a2*(x[k*ag+4*sg+r*rg+13]*alpha_i+x[k*ag+4*sg+r*rg+14])+a3*(x[k*ag+4*sg+r*rg+18]*alpha_i+x[k*ag+4*sg+r*rg+19]))*(qp*(cnt_matrix_cwc_p_v[a][k]+qc*cnt_matrix_cwc_c_v[a][k]))*N_tot_c_v_inv +(x[k*ag+5*sg+r*rg+3]*alpha_i+x[k*ag+5*sg+r*rg+4]+a1*(x[k*ag+5*sg+r*rg+8]*alpha_i+x[k*ag+5*sg+r*rg+9])+a2*(x[k*ag+5*sg+r*rg+13]*alpha_i+x[k*ag+5*sg+r*rg+14])+a3*(x[k*ag+5*sg+r*rg+18]*alpha_i+x[k*ag+5*sg+r*rg+19]))*(qp*(cnt_matrix_cwn_p_v[a][k]+qc*cnt_matrix_cwn_c_v[a][k]))*N_tot_n_v_inv;
                  
                  I_temp_n_v += (x[k*ag+3*sg+r*rg+3]*alpha_i+x[k*ag+3*sg+r*rg+4]+a1*(x[k*ag+3*sg+r*rg+8]*alpha_i+x[k*ag+3*sg+r*rg+9])+a2*(x[k*ag+3*sg+r*rg+13]*alpha_i+x[k*ag+3*sg+r*rg+14])+a3*(x[k*ag+3*sg+r*rg+18]*alpha_i+x[k*ag+3*sg+r*rg+19]))*(qp*(cnt_matrix_nwp_p_v[a][k]+qc*cnt_matrix_nwp_c_v[a][k]))*N_tot_p_v_inv + (x[k*ag+4*sg+r*rg+3]*alpha_i+x[k*ag+4*sg+r*rg+4]+a1*(x[k*ag+4*sg+r*rg+8]*alpha_i+x[k*ag+4*sg+r*rg+9])+a2*(x[k*ag+4*sg+r*rg+13]*alpha_i+x[k*ag+4*sg+r*rg+14])+a3*(x[k*ag+4*sg+r*rg+18]*alpha_i+x[k*ag+4*sg+r*rg+19]))*(qp*(cnt_matrix_nwc_p_v[a][k]+qc*cnt_matrix_nwc_c_v[a][k]))*N_tot_c_v_inv + (x[k*ag+5*sg+r*rg+3]*alpha_i+x[k*ag+5*sg+r*rg+4]+a1*(x[k*ag+5*sg+r*rg+8]*alpha_i+x[k*ag+5*sg+r*rg+9])+a2*(x[k*ag+5*sg+r*rg+13]*alpha_i+x[k*ag+5*sg+r*rg+14])+a3*(x[k*ag+5*sg+r*rg+18]*alpha_i+x[k*ag+5*sg+r*rg+19]))*(qp*(cnt_matrix_nwn_p_v[a][k]+qc*cnt_matrix_nwn_c_v[a][k]))*N_tot_n_v_inv;
              }
          }

          if (direct)
          {
              I_temp_p_v = I_temp_p; I_temp_c_v = I_temp_c; I_temp_n_v = I_temp_n;
          }
          
          int pj = max(ag*(a-1),0);
          int cj = max(ag*(a),0);
          double muBp, muBc, muBn, muBpv, muBcv, muBnv, mu_mat, cl, Icp, Icc, up, kpc, kpd, pro, pro_v_in, pro_v_out, rp;
          double p_vulp, p_vulc, p_vuln, p_vulpv, p_vulcv, p_vulnv;
          double xi_bp, xi_bc, xi_bn, xi_bpv, xi_bcv, xi_bnv;

          double ej1  = eta[a+1];
          double ej   = eta[a];
          double lossP, lossMS0, lossMS1;
          
          // Birth rate into each social group
       //   Rcpp::Rcout << "Birth rate into each social group: " << a << std::endl;

          if (a == 0)
          {
             // p_vulp = 0; p_vulc = sum_wcb; p_vuln = sum_wcb; p_vulpv = 0; p_vulcv = sum_wcb_v; p_vulnv = sum_wcb_v;
              p_vulp = 0; p_vulc = 0; p_vuln = 0; p_vulpv = 0; p_vulcv = 0; p_vulnv = 0;

              muBp = 0;    muBc = dailyBirthRate * (cov_c);   muBn = dailyBirthRate * (1 - cov_c);
              muBpv = 0;  muBcv = dailyBirthRate * (cov_c);  muBnv = dailyBirthRate * (1 - cov_c);
              xi_bp = 1;    xi_bc = 1;   xi_bn = 1;
              xi_bpv = 1;   xi_bcv = xi_boost;  xi_bnv = 1;
          }
          else
          {
              p_vulp = 0; p_vulc = 0; p_vuln = 0; p_vulpv = 0; p_vulcv = 0; p_vulnv = 0;
              muBp = 0; muBc = 0; muBn = 0;
              muBpv = 0; muBcv = 0; muBnv = 0;
              xi_bp = 1;    xi_bc = 1;   xi_bn = 1;
              xi_bpv = 1;   xi_bcv = xi_boost;  xi_bnv = 1;
          }
          

          for(int s = 0; s < 6; s++)
          {
            //  Rcpp::Rcout << "Social groups parameters: " << s << std::endl;
              
              kpd = kpc = pro = pro_v_in = pro_v_out = 0;
              cl = 0;     Icp = 1;    Icc = 0;    up = 0;
              
              double u, In, mu, vac_c_o, vac_c_i;
              int kp;
              if (a < 12)
              {
                  if (s == 0)
                  {   kp = 0*sg;       u = 0;             In = I_temp_p;   mu = muBp;  mu_mat = 0; p_vul = p_vulp; xi_b = xi_bp;
                  }
                  else if (s == 1)
                  {   kp = 1*sg;       u = (cov_c);       In = I_temp_c;   mu = muBc;  mu_mat = 0; p_vul = p_vulc; xi_b = xi_bc;
                  }
                  else if (s == 2)
                  {   kp = 2*sg;       u = (1-cov_c);     In = I_temp_n;   mu = muBn;  mu_mat = 0; p_vul = p_vuln; xi_b = xi_bn;
                  }
                  else if (s == 3)
                  {   kp = 3*sg;       u = 0;             In = I_temp_p_v; mu = muBpv; mu_mat = 0; p_vul = p_vulpv; xi_b = xi_bpv;
                  }
                  else if (s == 4)
                  {   kp = 4*sg;       u = (cov_c);       In = I_temp_c_v; mu = muBcv*(1-vac_cal(t1, 0)); mu_mat = muBcv*vac_cal(t1, 0); p_vul = p_vulcv;  xi_b = xi_bcv;
   
                  }
                  else
                  {   kp = 5*sg;       u = (1-cov_c);     In = I_temp_n_v; mu = muBnv; mu_mat = 0; p_vul = p_vulnv;  xi_b = xi_bnv;
                  }
              }
              else
              {
                  if (s == 0)
                  {   kp = 0*sg;       u = p_mat[a]*(1-cov_c);        In = I_temp_p;   mu = mu_mat = 0;  p_vul = p_vulp;
                  }
                  else if (s == 1)
                  {   kp = 1*sg;       u = p_mat[a]*(cov_c);          In = I_temp_c;   mu = mu_mat = 0;  p_vul = p_vulc;
                  }
                  else if (s == 2)
                  {   kp = 2*sg;       u = (1-p_mat[a]);              In = I_temp_n;   mu = mu_mat = 0;  p_vul = p_vuln;
                  }
                  else if (s == 3)
                  {   kp = 3*sg;       u = p_mat[a]*(1-cov_c);        In = I_temp_p_v; mu = mu_mat = 0; p_vul = p_vulpv;
                  }
                  else if (s == 4)
                  {   kp = 4*sg;       u = p_mat[a]*(cov_c);          In = I_temp_c_v; mu = mu_mat = 0; p_vul = p_vulcv;
                  }
                  else
                  {   kp = 5*sg;       u = (1-p_mat[a]);             In = I_temp_n_v; mu = mu_mat = 0; p_vul = p_vulnv;
                  }
              }
              
              for (int r = 0; r < 3; r++)
              {
               //   Rcpp::Rcout << "Risk groups parameters: " << r << std::endl;

                  double PST = 0;
                  x_tot = x_tot_1 = x_tot_2 = 0;
                  if (r==0){
                      rp = pVHR[a];
                  }
                  else if (r==1){
                      rp = pHR[a];}
                  else if (r==2){
                      rp = pLR[a];}
                  else
                      cout << "ERROR" << endl;
                  
                  if (r == 2)
                  {
                   //   Rcpp::Rcout << "Place where r == 2." << std::endl;
                      if (a==18 || a==19 || a==20)
                      {
                          if (a==18)
                              up = u18p;
                          else if (a==19)
                              up = u19p;
                          else if (a==20)
                              up = u20p;
                          else
                              cout << "ERROR" << endl;
                      }
                  }
              
                  if (s < 3)
                  {
                    //  Rcpp::Rcout << "Place where s < 3." << std::endl;
                      if (a < 12)
                      {
                          /*if (a > 0)
                              rp = 1;
                          
                          if (r == 0)
                              r_prop = VHR_g[a];
                          else if (r == 2)
                              r_prop = LR_g[a];
                          else
                              r_prop = 1;
                          
                          for (int i = 0; i < 24; i++)
                              PS[i] = (x[pj+0*sg+r*rg+i]+x[pj+1*sg+r*rg+i]+x[pj+2*sg+r*rg+i])*r_prop;*/
                          for (int i = 0; i < 21; i++)
                              PS[i] = (x[pj+0*sg+0*rg+i]+x[pj+0*sg+1*rg+i]+x[pj+0*sg+2*rg+i] + x[pj+1*sg+0*rg+i]+x[pj+1*sg+1*rg+i]+x[pj+1*sg+2*rg+i] + x[pj+2*sg+0*rg+i] + x[pj+2*sg+1*rg+i]+x[pj+2*sg+2*rg+i]);
                          
                          for (int i = 21; i < 24; i++)
                              PS[i] = x[pj+0*sg+r*rg+i] + x[pj+1*sg+r*rg+i] + x[pj+2*sg+r*rg+i];
                      }
                      else
                      {
                          for (int i = 0; i < 24; i++)
                              PS[i] = (x[pj+0*sg+0*rg+i]+x[pj+0*sg+1*rg+i]+x[pj+0*sg+2*rg+i] + x[pj+1*sg+0*rg+i]+x[pj+1*sg+1*rg+i]+x[pj+1*sg+2*rg+i] + x[pj+2*sg+0*rg+i] + x[pj+2*sg+1*rg+i]+x[pj+2*sg+2*rg+i]);
                      }
                  }
                  else
                  {
                   //   Rcpp::Rcout << "Place where s >= 3." << std::endl;

                      if (a < 12)
                      {
                         // Rcpp::Rcout << "Place where a < 12." << std::endl;
                          /*if (a > 0)
                              rp = 1;
                          
                          if (r == 3)
                              r_prop = VHR_g[a];
                          else if (r == 5)
                              r_prop = LR_g[a];
                          else
                              r_prop = 1;
                          
                          for (int i = 0; i < 24; i++)
                              PS[i] = (x[pj+3*sg+r*rg+i]+x[pj+4*sg+r*rg+i]+x[pj+5*sg+r*rg+i])*r_prop;
                           */
                          for (int i = 0; i < 21; i++)
                              PS[i] = x[pj+3*sg+0*rg+i] + x[pj+3*sg+1*rg+i] + x[pj+3*sg+2*rg+i] + x[pj+4*sg+0*rg+i] + x[pj+4*sg+1*rg+i] + x[pj+4*sg+2*rg+i] + x[pj+5*sg+0*rg+i] + x[pj+5*sg+1*rg+i] + x[pj+5*sg+2*rg+i];
                          
                          for (int i = 21; i < 24; i++)
                              PS[i] = x[pj+3*sg+r*rg+i] + x[pj+4*sg+r*rg+i] + x[pj+5*sg+r*rg+i];
 
                      }
                      else
                      {
                     //     Rcpp::Rcout << "Place where a >= 12." << std::endl;

                          for (int i = 0; i < 24; i++)
                              PS[i] = x[pj+3*sg+0*rg+i] + x[pj+3*sg+1*rg+i] + x[pj+3*sg+2*rg+i] + x[pj+4*sg+0*rg+i] + x[pj+4*sg+1*rg+i] + x[pj+4*sg+2*rg+i] + x[pj+5*sg+0*rg+i] + x[pj+5*sg+1*rg+i] + x[pj+5*sg+2*rg+i];
                      }
                  }
                  
                  double cpmu, cpo, cMmu, cMo, cLo, cMa, cpmu_dose, cpo_dose, cMmu_dose, cMo_dose, cLo_dose, cMa_dose;

                  int p = a*ag + s*sg + r*rg;
                  int q = a*ag + kpc + r*rg;
                  int o = 0;
                  if (s>3)
                      o = a*ag + (s-3)*sg + r*rg;
                  
                  if (s < 3)
                  {

                      cpmu = 0; cpmu_dose = 0;
                      cpo = 0; cpo_dose = 0;
                      cMmu = 0; cMmu_dose = 0;
                      cMo = 0; cMo_dose = 0;
                      cLo = 0; cLo_dose = 0;
                      cMa = 0; cMa_dose = 0;
                      SN = 1;
                      lossMS1 = lossMS0 = lossP = 0;

                  }
                  else
                  {

                      lossP = x[a*ag + s*sg + r*rg + 21]*(1.0/60.0);
                      lossMS0 = x[a*ag + s*sg + r*rg + 22]*(om_mab);
                      lossMS1 = x[a*ag + s*sg + r*rg + 23]*(om_mab);
                      
                      cpmu = 0; cpmu_dose = 0;
                      cpo = 0; cpo_dose = 0;
                      cMmu = 0; cMmu_dose = 0;
                      cMo = 0; cMo_dose = 0;
                      cLo = 0; cLo_dose = 0;
                      cMa = 0; cMa_dose = 0;

                      if (a == 0){
                      //    Rcpp::Rcout << "Place where s >= 3 and a == 0." << std::endl;

                          if (r == 0){
                     //         Rcpp::Rcout << "Place where s >= 3 and a == 0 and r == 0." << std::endl;
                              cpo = cal_pal(t1, 0);
                              cMo = cal_mAB_VHR(t1, 0);
                              cpo_dose = cal_pal_dose(t1, 0);
                              cMo_dose = cal_mAB_VHR_dose(t1, 0);}
                          else if (r == 1){
                         //     Rcpp::Rcout << "Place where s >= 3 and a == 0 and r == 1." << std::endl;
                              cMo = cal_mAB_HR(t1, 0);
                              cMo_dose = cal_mAB_HR_dose(t1, 0);
                          }
                          else{
                          //    Rcpp::Rcout << "Place where s >= 3 and a == 0 and r == 2." << std::endl;
                              cMo = cal_mAB_LR(t1, 0);
                              cMo_dose = cal_mAB_LR_dose(t1, 0);
                          }
                      }
                      else if (a > 0)
                      {
                          //Rcpp::Rcout << "Place where s >= 3 and a > 0." << std::endl;

                          if (r == 0){
                           //   Rcpp::Rcout << "Place where s >= 3 and a > 0 and r == 0." << std::endl;
                              cpo = cal_pal(t1, a);
                              cpo_dose = cal_pal_dose(t1, a);
                              cMo = cal_mAB_VHR(t1, a);
                              cMo_dose = cal_mAB_VHR_dose(t1, a);
                              S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];
                              SN = 1;
                          }
                          else if (r == 1){
                          //    Rcpp::Rcout << "Place where s >= 3 and a > 0 and r == 1." << std::endl;
                              cMo = cal_mAB_HR(t1, a); // Rcpp::Rcout << "A" << endl;
                              cLo = cal_LAV_HR(t1, a);// Rcpp::Rcout << "B" << endl;
                              cMo_dose = cal_mAB_HR_dose(t1, a);// Rcpp::Rcout << "C" << endl;
                              cLo_dose = cal_LAV_HR_dose(t1, a);// Rcpp::Rcout << "D" << endl;
                              S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16]; //Rcpp::Rcout << "E" << endl;
                              SN = 1;// Rcpp::Rcout << "F" << endl;
                          //    Rcpp::Rcout << "Place where s >= 3 and a > 0 and r == 1. END" << std::endl;
                          }
                          else
                          {
                          //    Rcpp::Rcout << "Place where s >= 3 and a > 0 and r == 2." << std::endl;
                              cMo = cal_mAB_LR(t1, a);
                              cMo_dose = cal_mAB_LR_dose(t1, a);
                              if (s==4)
                              {
                           //       Rcpp::Rcout << "Place where s == 4 and a > 0 and r == 2." << std::endl;
                                  cMa = vac_cal(t1, a);
                                  cMa_dose = vac_cal_dose(t1, a);
                                  S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];
                              }
                              else
                              {
                           //       Rcpp::Rcout << "Place where s == 3 or 5 and a > 0 and r == 2." << std::endl;
                                  cMa = 0;
                                  cMa_dose = 0;
                                  SN = 1;
                                  cLo = cal_LAV_LR(t1, a);
                                  cLo_dose = cal_LAV_LR_dose(t1, a);
                                  S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];
                              }
                          }
                      }
                      else
                      {
                          
                      }
                  }
                  
                  for (int i = 0; i < 21; i++)
                  {
                      x_tot += x[o+i];
                      PST += PS[i];
                      if (PST < 1)
                          PST = 1;
                  }
                  for (int i = 0; i < 2; i++)
                      x_tot_1 += x[o+i];
                  
                  for (int i = 2; i < 6; i++)
                      x_tot_2 += x[o+i];
                  
              //    Rcpp::Rcout << "keeping track of doses below" << std::endl;

                  protectpal += cpo_dose*x_tot;
                  protectmabs += cMo_dose*x_tot;
                  protectLAV += cLo_dose*x_tot;
                  protectmat += cMa_dose*x_tot;

                  dxdt[p+0] = (1.0-p_vul)*mu*rp + mu_mat*rp - x[p+0]*xi*xi_b - (x[p+0])*ej1 + PS[0]*ej*rp*u - x[o+0]*(cMo) - x[o+0]*(cpo);
                  
                  dxdt[p+1] = p_vul*mu*rp  + x[p+0]*xi*xi_b + lossP + lossMS0 - x[p+1]*In*beta - (x[p+1])*ej1 + PS[1]*ej*rp*u - x[o+1]*(cMo) - x[o+1]*(cpo) - cLo*S0 - cMa*S0;
                  dxdt[p+2] = x[p+1]*In*beta                 - x[p+2]*si          - (x[p+2])*ej1 + PS[2]*ej*rp*u - x[o+2]*(cMo) - x[o+2]*(cpo);
                  dxdt[p+3] = x[p+2]*si*pA[a]                - x[p+3]*ga0*rho     - (x[p+3])*ej1 + PS[3]*ej*rp*u - x[o+3]*(cMo) - x[o+3]*(cpo);
                  dxdt[p+4] = x[p+2]*si*(1.0-pA[a])          - x[p+4]*ga0         - (x[p+4])*ej1 + PS[4]*ej*rp*u - x[o+4]*(cMo) - x[o+4]*(cpo);
                  dxdt[p+5] = x[p+4]*ga0 + x[p+3]*ga0*rho    - x[p+5]*om          - (x[p+5])*ej1 + PS[5]*ej*rp*u - x[o+5]*(cMo) - x[o+5]*(cpo) + cLo*S0 + cMa*S0;
                  
                  dxdt[p+6] = x[p+5]*om                      - d1*x[p+6]*In*beta + lossMS1  - (x[p+6])*ej1 + PS[6]*ej*rp*u - x[o+6]*(cMo) - x[o+6]*(cpo) - cLo*S1 - cMa*S1;
                  dxdt[p+7] = d1*x[p+6]*In*beta              - x[p+7]*si          - (x[p+7])*ej1 + PS[7]*ej*rp*u - x[o+7]*0 - x[o+7]*(cpo);
                  dxdt[p+8] = x[p+7]*si*pA[a]                - x[p+8]*ga1*rho     - (x[p+8])*ej1 + PS[8]*ej*rp*u - x[o+8]*0 - x[o+8]*(cpo);
                  dxdt[p+9] = x[p+7]*si*(1.0-pA[a])          - x[p+9]*ga1         - (x[p+9])*ej1 + PS[9]*ej*rp*u - x[o+9]*0 - x[o+9]*(cpo);
                  dxdt[p+10] = x[p+9]*ga1 + x[p+8]*ga1*rho   - x[p+10]*om         - (x[p+10])*ej1 + PS[10]*ej*rp*u - x[o+10]*0 - x[o+10]*(cpo) + cLo*S1 + cMa*S1;
                  
                  dxdt[p+11] = x[p+10]*om                    - d2*x[p+11]*In*beta - (x[p+11])*ej1 + PS[11]*ej*rp*u - x[o+11]*0 - x[o+11]*(cpo) - cLo*S2 - cMa*S2;
                  dxdt[p+12] = d2*x[p+11]*In*beta            - x[p+12]*si         - (x[p+12])*ej1 + PS[12]*ej*rp*u - x[o+12]*0 - x[o+12]*(cpo);
                  dxdt[p+13] = x[p+12]*si*pA[a]              - x[p+13]*ga2*rho    - (x[p+13])*ej1 + PS[13]*ej*rp*u - x[o+13]*0 - x[o+13]*(cpo);
                  dxdt[p+14] = x[p+12]*si*(1.0-pA[a])        - x[p+14]*ga2        - (x[p+14])*ej1 + PS[14]*ej*rp*u - x[o+14]*0 - x[o+14]*(cpo);
                  dxdt[p+15] = x[p+14]*ga2 + x[p+13]*ga2*rho - x[p+15]*om         - (x[p+15])*ej1 + PS[15]*ej*rp*u - x[o+15]*0 - x[o+15]*(cpo) + cLo*S2 + cMa*S2;
                  
                  dxdt[p+16] = x[p+15]*om + x[p+20]*om       - d3*x[p+16]*In*beta - (x[p+16])*ej1 + PS[16]*ej*rp*u - x[o+16]*0 - x[o+16]*(cpo) - cLo*S3 - cMa*S3;
                  dxdt[p+17] = d3*x[p+16]*In*beta            - x[p+17]*si         - (x[p+17])*ej1 + PS[17]*ej*rp*u - x[o+17]*0 - x[o+17]*(cpo);
                  dxdt[p+18] = x[p+17]*si*pA[a]              - x[p+18]*ga3*rho    - (x[p+18])*ej1 + PS[18]*ej*rp*u - x[o+18]*0 - x[o+18]*(cpo);
                  dxdt[p+19] = x[p+17]*si*(1.0-pA[a])        - x[p+19]*ga3        - (x[p+19])*ej1 + PS[19]*ej*rp*u - x[o+19]*0 - x[o+19]*(cpo);
                  dxdt[p+20] = x[p+19]*ga3 + x[p+18]*ga3*rho - x[p+20]*om         - (x[p+20])*ej1 + PS[20]*ej*rp*u - x[o+20]*0 - x[o+20]*(cpo) + cLo*S3 + cMa*S3;
                  
                  // Vaccination groups
                  dxdt[p+21] = x_tot*cpo - lossP - x[p+21]*ej1 + PS[21]*ej*u*rp;  // pal or mabs
                  dxdt[p+22] = x_tot_1*cMo - lossMS0 - x[p+22]*ej1 + PS[22]*ej*u*rp; // + x[(a-1)*ag + s*sg + r*rg +22]*ej;
                  dxdt[p+23] = x_tot_2*cMo - lossMS1 - x[p+23]*ej1 + PS[23]*ej*u*rp; // + x[(a-1)*ag +

                 // dxdt[p+21] = mu*rp*cpmu + cpo*PST     - x[p+21]*ej1 + x[p+21]*ej1;
                 // dxdt[p+22] = mu*rp*cMmu + PST2*ej*u*rp*cMo - loss    - x[p+22]*ej1 + PS[22]*ej*u*rp*(1-cpo)*(1-cMo) ; // + x[(a-1)*ag + s*sg + r*rg +22]*ej;
              }
          }
          // Vaccine groups
          dxdt[a*ag + 6*sg + 0] = si*(x[cj+3*sg+0*rg+2] + x[cj+4*sg+0*rg+2] + x[cj+5*sg+0*rg+2] + x[cj+3*sg+1*rg+2] + x[cj+4*sg+1*rg+2] + x[cj+5*sg+1*rg+2] + x[cj+3*sg+2*rg+2] + x[cj+4*sg+2*rg+2] + x[cj+5*sg+2*rg+2]);
          dxdt[a*ag + 6*sg + 1] = protectpal;
          dxdt[a*ag + 6*sg + 2] = protectmabs;
         // cout << protectmabs << endl;
          /*dxdt[a*ag + 6*sg + 2] = protectmabs + vac_cal(t1, a)*(x[a*ag+4*sg+2*rg+1]+x[a*ag+4*sg+2*rg+6]+x[a*ag+4*sg+2*rg+11]+x[a*ag+4*sg+2*rg+16])
          + (cal_LAV_HR(t1, a))*((x[a*ag+0*sg+1*rg+1]+x[a*ag+0*sg+1*rg+6]+x[a*ag+0*sg+1*rg+11]+x[a*ag+0*sg+1*rg+16]) + (x[a*ag+1*sg+1*rg+1]+x[a*ag+1*sg+1*rg+6]+x[a*ag+1*sg+1*rg+11]+x[a*ag+1*sg+1*rg+16]))
          + (cal_LAV_LR(t1, a))*((x[a*ag+0*sg+2*rg+1]+x[a*ag+0*sg+2*rg+6]+x[a*ag+0*sg+2*rg+11]+x[a*ag+0*sg+2*rg+16]) + (x[a*ag+1*sg+2*rg+1]+x[a*ag+1*sg+2*rg+6]+x[a*ag+1*sg+2*rg+11]+x[a*ag+1*sg+2*rg+16]));
          */
          dxdt[a*ag + 6*sg + 3] = protectLAV;
          
          dxdt[a*ag + 6*sg + 4] = protectmat;
          dxdt[a*ag + 6*sg + 5] = 0;
          dxdt[a*ag + 6*sg + 6] = 0;
          dxdt[a*ag + 6*sg + 7] = 0;
          
          // Monitoring Parent
          dxdt[a*ag + 6*sg + 8] =  si*(x[cj+3*sg+0*rg+2]+x[cj+3*sg+0*rg+7]+x[cj+3*sg+0*rg+12]+x[cj+3*sg+0*rg+17]); //VHR
          dxdt[a*ag + 6*sg + 9] =  si*(x[cj+3*sg+1*rg+2]+x[cj+3*sg+1*rg+7]+x[cj+3*sg+1*rg+12]+x[cj+3*sg+1*rg+17]); //HR
          dxdt[a*ag + 6*sg + 10] = si*(x[cj+3*sg+2*rg+2]+x[cj+3*sg+2*rg+7]+x[cj+3*sg+2*rg+12]+x[cj+3*sg+2*rg+17]); //LR
          
          // Monitoring Cocoon
          dxdt[a*ag + 6*sg + 11] = si*(x[cj+4*sg+0*rg+2]+x[cj+4*sg+0*rg+7]+x[cj+4*sg+0*rg+12]+x[cj+4*sg+0*rg+17]); //VHR
          dxdt[a*ag + 6*sg + 12] = si*(x[cj+4*sg+1*rg+2]+x[cj+4*sg+1*rg+7]+x[cj+4*sg+1*rg+12]+x[cj+4*sg+1*rg+17]); //HR
          dxdt[a*ag + 6*sg + 13] = si*(x[cj+4*sg+2*rg+2]+x[cj+4*sg+2*rg+7]+x[cj+4*sg+2*rg+12]+x[cj+4*sg+2*rg+17]); //LR
          
          // Monitoring Neither
          dxdt[a*ag + 6*sg + 14] = si*(x[cj+5*sg+0*rg+2]+x[cj+5*sg+0*rg+7]+x[cj+5*sg+0*rg+12]+x[cj+5*sg+0*rg+17]); //VHR
          dxdt[a*ag + 6*sg + 15] = si*(x[cj+5*sg+1*rg+2]+x[cj+5*sg+1*rg+7]+x[cj+5*sg+1*rg+12]+x[cj+5*sg+1*rg+17]); //HR
          dxdt[a*ag + 6*sg + 16] = si*(x[cj+5*sg+2*rg+2]+x[cj+5*sg+2*rg+7]+x[cj+5*sg+2*rg+12]+x[cj+5*sg+2*rg+17]); //LR
          
          x[a*ag + 6*sg + 17] = N_tot_n[a];
          x[a*ag + 6*sg + 18] = N_tot_c[a];
          x[a*ag + 6*sg + 19] = N_tot_p[a];
          
          x[a*ag + 6*sg + 20] = N_tot_n_v[a];
          x[a*ag + 6*sg + 21] = N_tot_c_v[a];
          x[a*ag + 6*sg + 22] = N_tot_p_v[a];
      }
  }
};



List RunInterventions::SampleWeekly(List vac_calendar, List vac_dose, double cov_c, List vac_info, VectorXd posteriors)
{
  // Restart values
  this->currentODETime = this->run_start;
  VectorXd inc_tot(this->A);
  inc_tot = VectorXd::Zero(this->A);
  this->weekNo = 0;
  this->dayNoAfterBurn = 0;
    
  VectorXd currentParamValues = posteriors;
  ParameterValuesforODE(currentParamValues);

  // Set up and Run ODE solver
  EulerT<state_t> integrator;
  SystemT<state_t, system_t> System;
  asc::Recorder recorder;
  vector< double >  x0 = generateInitialStates(cov_c);
  ODE_desc ODE_desc_inst(this, vac_calendar, vac_dose, vac_info, cov_c);
// 573
    NumericMatrix sampleWeeklyIncidence(521, this->A*9);
    NumericMatrix no_doses(521, 4 * 25);

    while (this->currentODETime < (this->run_full + this->run_burn)){

      integrator(ODE_desc_inst, x0, this->currentODETime, this->dt);
      if (this->currentODETime > this->run_burn){

        getWeeklyIncidence(x0, sampleWeeklyIncidence, no_doses, FALSE);
          
      }
    }
    return Rcpp::List::create(Rcpp::Named("inci") = sampleWeeklyIncidence,
                              Rcpp::Named("doses") = no_doses);
}


List RunInterventions::SampleMonthly(List vac_calendar, List vac_dose, double cov_c, List vac_info, VectorXd posteriors)
{
  // Restart values
  this->currentODETime = this->run_start;
  VectorXd inc_tot(this->A);
  inc_tot = VectorXd::Zero(this->A);
  this->monthNo = 0;
  this->dayNoAfterBurn = 0;
    
  VectorXd currentParamValues = posteriors;
  ParameterValuesforODE(currentParamValues);

  // Set up and Run ODE solver
  EulerT<state_t> integrator;
  SystemT<state_t, system_t> System;
  asc::Recorder recorder;
  vector< double >  x0 = generateInitialStates(cov_c);
  ODE_desc ODE_desc_inst(this, vac_calendar, vac_dose, vac_info, cov_c);

    NumericMatrix sampleMonthlyIncidence(120, this->A*9);
    NumericMatrix no_doses(120, 4 * this->A);

    while (this->currentODETime < (this->run_full + this->run_burn)){

      integrator(ODE_desc_inst, x0, this->currentODETime, this->dt);
      if (this->currentODETime > (this->run_burn)){
        getMonthlyIncidence(x0, sampleMonthlyIncidence, no_doses, FALSE);
          
      }
    }
    return Rcpp::List::create(Rcpp::Named("inci") = sampleMonthlyIncidence,
                              Rcpp::Named("doses") = no_doses);
}

NumericMatrix RunInterventions::getDebug(List vac_calendar, List vac_dose, double cov_c, List vac_info, VectorXd currentParamValues)
{
  // Restart values
  this->currentODETime = 0;
  this->dayNoAfterBurn = 0;
  ParameterValuesforODE(currentParamValues);

  // Set up and Run ODE solver
  EulerT<state_t> integrator;
  SystemT<state_t, system_t> System;
  asc::Recorder recorder;
  vector< double > x0 = generateInitialStates(cov_c);
  ODE_desc ODE_desc_inst(this, vac_calendar, vac_dose, vac_info, cov_c);
  
  NumericMatrix sampleDailyIncidence(365, 11375);
  
  while (this->currentODETime < 365){
     for (int k = 0; k < 11375; k++) {
        sampleDailyIncidence(this->currentODETime, k) = x0[k];
      }
      integrator(ODE_desc_inst, x0, this->currentODETime, this->dt);

    }
  return sampleDailyIncidence; // RETURN THE LOG LIKELIHOOD
}

#endif /* RunIntervention_h */
