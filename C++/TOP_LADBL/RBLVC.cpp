#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List RBLVC (arma::vec x, arma::vec y, arma:: mat w, arma:: mat c, arma::mat e, int maxSteps,int n, arma::vec hatAlpha,arma::vec hatb, double hatBeta, arma::vec hatEta, double hatTau, arma::vec hatV, double hatSg1, arma:: vec hatSg2, arma::mat invSigAlpha0, arma:: mat invSigb0, double hatEtaSq1,double hatEtaSq2, double xi1, double xi2, double r1, double r, double a, double b, int progress)
{
    unsigned int q1 = e.n_cols, q2 = c.n_cols;
    arma::mat gsAlpha(maxSteps, q1),
            gsb(maxSteps,q2),
            gsEta(maxSteps,q1),
            gsV(maxSteps, n),
            gsSg2(maxSteps, q1);
    
    arma::vec gsBeta(maxSteps),
              gsEtaSq1(maxSteps),
              gsEtaSq2(maxSteps),
              gsSg1(maxSteps),
              gsTau(maxSteps);
           
    
    arma::mat varAlpha,varb, tEEoV(q1,q1),tCCoV(q2,q2);
    arma::vec res, REoV(q1),RCoV(q2), meanAlpha, muV, muS2, meanb;

    double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2), XgXgoV1,XgXgoV2, RXgoV1, RXgoV2, meanG1, varG1,meanG2,varG2, ResSqoV, muS1;
    
    
    for (int k = 0; k < maxSteps; k++) {
        // Rcpp::Rcout << "alpha" << std::endl;
        res = y - c*hatb - x * hatBeta - xi1*hatV -w*hatEta;
        tEEoV = (e.each_col()/hatV).t() * e;
        REoV = arma::sum(e.each_col()% (res/hatV), 0).t();
        varAlpha = arma::inv_sympd(tEEoV*hatTau/xi2Sq + invSigAlpha0);
        meanAlpha = varAlpha * REoV * hatTau / xi2Sq;
        hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
        res -= e * hatAlpha;
        gsAlpha.row(k) = hatAlpha.t();
        
        
        // Rcpp::Rcout << "alpha" << std::endl;
        
        tCCoV = (c.each_col()/hatV).t() * c;
        res += c*hatb;
        RCoV = arma::sum(c.each_col()% (res/hatV), 0).t();
        varb = arma::inv_sympd(tCCoV*hatTau/xi2Sq + invSigb0);
        meanb = varb * RCoV * hatTau / xi2Sq;
        hatb = mvrnormCpp(meanb, varb);
        res -= c* hatb;
        gsb.row(k) = hatb.t();
        
        // Rcpp::Rcout << "v" << std::endl;
        res += xi1*hatV;
        lambV = hatTau*xi1Sq/xi2Sq + 2*hatTau;
        muV = arma::sqrt((xi1Sq+2*xi2Sq) / arma::square(res));
        for(unsigned int i = 0; i<n; i++){
            bool flag = true;
            while(flag){
                hatV(i) = 1/rinvGauss(muV(i), lambV);
                if(hatV(i)<=0 || std::isinf(hatV(i)) || std::isnan(hatV(i))){
                    if(progress != 0) Rcpp::Rcout << "hatV(i) <= 0 or nan or inf" << std::endl;
                    Rcpp::checkUserInterrupt();
                }else{
                    flag = false;
                }
            }
        }
        res -= xi1*hatV;
        gsV.row(k) = hatV.t();
        
        
        // Rcpp::Rcout << "S1" << std::endl;
        muS1 = std::sqrt(hatEtaSq1/ pow(hatBeta,2));

        bool flag = true;
        while(flag){
            hatSg1 = 1/rinvGauss(muS1, hatEtaSq1);
            if(hatSg1<=0 || std::isinf(hatSg1) || std::isnan(hatSg1)){
                if(progress != 0) Rcpp::Rcout << "hatSg1: " << hatSg1 << std::endl;
                    Rcpp::checkUserInterrupt();
                }else{
                    flag = false;
                }
        }
        
        gsSg1(k) = hatSg1;
        
        res += x * hatBeta;
        XgXgoV1 = arma::as_scalar((x/hatV).t() * x);
        varG1 = 1/(XgXgoV1*hatTau/xi2Sq + 1/hatSg1);
                   
        RXgoV1 = arma::sum(x % (res/hatV));
        meanG1 = varG1 * RXgoV1 * hatTau / xi2Sq;
        hatBeta = R::rnorm(meanG1, sqrt(varG1));
        res -= x * hatBeta;
        
        gsBeta(k) = hatBeta;
    
    
       muS2 = std::sqrt(hatEtaSq2)/ arma::abs(hatEta);
       for(unsigned int j = 0; j<q1; j++){
           bool flag = true;
           while(flag){
               hatSg2(j) = 1/rinvGauss(muS2(j), hatEtaSq2);
               if(hatSg2(j)<=0 || std::isinf(hatSg2(j)) || std::isnan(hatSg2(j))){
                   if(progress != 0) Rcpp::Rcout << "hatSg2(j): " << hatSg2(j) << std::endl;
                   Rcpp::checkUserInterrupt();
               }else{
                   flag = false;
               }
           }
       }
       gsSg2.row(k) = hatSg2.t();
        
        for(unsigned int j=0; j<q1; j++){
            res += w.col(j) * hatEta(j);
            XgXgoV2 = arma::as_scalar((w.col(j)/hatV).t() * w.col(j));
            varG2 = 1/(XgXgoV2*hatTau/xi2Sq + 1/hatSg2(j));
            
            RXgoV2 = arma::sum(w.col(j) % (res/hatV));
            meanG2 = varG2 * RXgoV2 * hatTau / xi2Sq;
            hatEta(j) = R::rnorm(meanG2, sqrt(varG2));
            res -= w.col(j) * hatEta(j);
        }
        gsEta.row(k) = hatEta.t();
    
    // Rcpp::Rcout << "etaSq1" << std::endl;
        double shape1 = 1+1;
        double rate1 = hatSg1/2 + r1;
        hatEtaSq1 = R::rgamma(shape1, 1/rate1);
        gsEtaSq1(k) = hatEtaSq1;
        
        
        // Rcpp::Rcout << "tau" << std::endl;
        double shape = a + 3*n/2;
        ResSqoV = arma::accu(arma::square(res)/hatV);
        double rate = b + arma::accu(hatV) + ResSqoV/(2*xi2Sq);
        hatTau = R::rgamma(shape, 1/rate);
        gsTau(k) = hatTau;
        
        // Rcpp::Rcout << "eta2Sq2" << std::endl;
        double shape2 = q1+1;
        double rate2 = arma::accu(hatSg2)/2 + r;
        hatEtaSq2 = R::rgamma(shape2, 1/rate2);
        gsEtaSq2(k) = hatEtaSq2;
        
        
    }
    
    return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,
                            Rcpp::Named("GS.b")=gsb,
                            Rcpp::Named("GS.beta") = gsBeta,
                            Rcpp::Named("GS.eta") = gsEta,
                            Rcpp::Named("GS.tau") = gsTau,
                            Rcpp::Named("GS.v") = gsV,
                            Rcpp::Named("GS.s1") = gsSg1,
                            Rcpp::Named("GS.s2") = gsSg2,
                            Rcpp::Named("GS.eta21.sq") = gsEtaSq1,
                            Rcpp::Named("GS.eta22.sq") = gsEtaSq2);
}
