#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BLSSVC (arma::vec x, arma::vec y, arma:: mat e, arma:: mat c, arma:: mat w, int maxSteps, int n, double hatBeta, arma:: vec hatEta, arma::vec hatAlpha,arma::vec hatb, double hatInvTauSq1, arma:: vec hatInvTauSq2,double sg1, arma::vec sg2,double hatPiEta, double hatPiBeta, arma::mat invSigAlpha0, arma:: mat invSigb0,double hatLambdaSqStar1, double hatLambdaSqStar2, double hatSigmaSq, double aStar, double bStar, double alpha, double gamma, double mu0, double nu0,int progress)
{
    unsigned int q1 = e.n_cols, q2 = c.n_cols;
    arma::mat gsAlpha(maxSteps, q1),
            gsb(maxSteps, q2),
            gseta(maxSteps,q1),
            gsInvTauSq2(maxSteps, q1),
            gsSS2(maxSteps,q1);
        
    arma::vec  gsBeta(maxSteps),
            gsInvTauSq1(maxSteps),
            gsLambdaStar1(maxSteps),
            gsLambdaStar2(maxSteps),
            gsSigmaSq(maxSteps),
            gsPiBeta(maxSteps),
            gsPiEta(maxSteps),
            gsSS1(maxSteps);


    arma::mat tEE= e.t()*e, tCC = c.t()*c, varAlpha, varb;
    arma::vec res, meanAlpha, meanb, tRsRs2, muInvTauSq2;
    double txx, tempS1, tempS2, varRs1, meanRs1, meanRs2, varRs2, lInvTauSq1, lB, lE,t, lInvTauSq2,muInvTauSq1,tRsRs1,GjtRes, WjtRes;
    
    arma:: mat tww = w.t()*w;
    arma::vec tBrBr2Diag = tww.diag();
    
    
    for (int k = 0; k < maxSteps; k++) {
        // alpha|
        varAlpha = arma::inv(tEE/hatSigmaSq + invSigAlpha0);
        res = y - c*hatb- x * hatBeta-w*hatEta;
        meanAlpha = varAlpha * (e.t() * res/hatSigmaSq);
        hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
        res -= e* hatAlpha;
        gsAlpha.row(k) = hatAlpha.t();
        
        // b|
        varb = arma::inv(tCC/hatSigmaSq + invSigb0);
        res += c*hatb;
        meanb = varb * (c.t() * res/hatSigmaSq);
        hatb = mvrnormCpp(meanb, varb);
        res -= c* hatb;
        gsb.row(k) = hatb.t();
        
        // Beta|
        
        txx = arma::accu(arma::square(x));
        tempS1 = 1/(txx + hatInvTauSq1);
        varRs1 = hatSigmaSq * tempS1;
        res += x * hatBeta;
        GjtRes = arma::as_scalar(x.t() * res);
        meanRs1 = tempS1 * GjtRes;
        lB = hatPiBeta/(hatPiBeta + (1-hatPiBeta)*sqrt(hatInvTauSq1*tempS1)*exp(0.5/hatSigmaSq*tempS1*pow(GjtRes,2)));
        
        t = R::runif(0, 1);
        if(t<lB){
            hatBeta = 0; sg1=0;
        }else{
            hatBeta = R::rnorm(meanRs1, sqrt(varRs1));sg1=1;
        }
        res -= x * hatBeta;
        gsBeta(k) = hatBeta;
        gsSS1(k) = sg1;
        
        // eta|
        
        for(unsigned int j=0;j<q1;j++){
            tempS2 = 1/(tBrBr2Diag(j) + hatInvTauSq2(j));
            varRs2 = hatSigmaSq * tempS2;
            res += w.col(j) * hatEta(j);
            WjtRes = arma::as_scalar(w.col(j).t() * res);
            meanRs2 = tempS2*WjtRes;
            
            lE = hatPiEta/(hatPiEta + (1-hatPiEta)*sqrt(hatInvTauSq2(j)*tempS2)*exp(0.5/hatSigmaSq*tempS2*pow(WjtRes,2)));
            
            t = R::runif(0, 1);
            if(t<lE){
                hatEta(j) = 0;sg2(j)=0;
            }else{
                hatEta(j) = R::rnorm(meanRs2, sqrt(varRs2));sg2(j)=1;
                                 
            }
            res -= w.col(j) * hatEta(j);
        }

        gseta.row(k) = hatEta.t();
        gsSS2.row(k) = sg2.t();
        
        // sigma.sq|
        double shapeSig = alpha + n/2 + sg1/2+ arma::accu(hatEta!=0)/2;
        double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) +
                                      pow(hatBeta,2)*hatInvTauSq1+ arma::accu(square(hatEta) % hatInvTauSq2));
        hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
        gsSigmaSq(k) = hatSigmaSq;
        
        
        // invTAUsq.star1|

        lInvTauSq1 = hatLambdaSqStar1;
        tRsRs1 = pow(hatBeta,2);
        muInvTauSq1 = sqrt(hatLambdaSqStar1 * hatSigmaSq / tRsRs1);
        if(hatBeta==0){
            hatInvTauSq1 = 1/R::rgamma(1,2/lInvTauSq1);
        }else{
            hatInvTauSq1 = rinvgaussian(muInvTauSq1, lInvTauSq1);
            
        }
        
        gsInvTauSq1(k) = hatInvTauSq1;
        
        // invTAUsq.star2|
    
        lInvTauSq2 = hatLambdaSqStar2;
        tRsRs2 = arma::square(hatEta);
        muInvTauSq2 = arma::sqrt(hatLambdaSqStar2 * hatSigmaSq / tRsRs2);
        for(unsigned int j = 0; j<q1; j++){
            if(hatEta(j)==0){
                hatInvTauSq2(j) = 1/R::rgamma(1,2/lInvTauSq2);
            }else{
                hatInvTauSq2(j) = rinvgaussian(muInvTauSq2(j), lInvTauSq2);
            }
            
        }
            
        gsInvTauSq2.row(k) = hatInvTauSq2.t();
        
        // lambda.star1|
        double shapeS1 = aStar + 1;
        double rateS1 = bStar + (1/hatInvTauSq1)/2;
        hatLambdaSqStar1 = R::rgamma(shapeS1, 1/rateS1);
        gsLambdaStar1(k) = hatLambdaSqStar1;
        
        // lambda.star2|
        double shapeS2 = aStar + q1;
        double rateS2 = bStar + arma::accu(1/hatInvTauSq2)/2;
        hatLambdaSqStar2 = R::rgamma(shapeS2, 1/rateS2);
        gsLambdaStar2(k) = hatLambdaSqStar2;
        
        double shape1_b = mu0 + 1-sg1;
        double shape2_b = nu0 + sg1;
        hatPiBeta = R::rbeta(shape1_b, shape2_b);
        gsPiBeta(k) = hatPiBeta;
        
        // pi.star|
        double shape1_e = mu0 + arma::accu(hatEta == 0);
        double shape2_e = nu0 + arma::accu(hatEta != 0);
        hatPiEta = R::rbeta(shape1_e, shape2_e);
        gsPiEta(k) = hatPiEta;
        
    }
    
    
    return Rcpp::List::create(Rcpp::Named("GS.b") = gsb,
                            Rcpp::Named("GS.alpha") = gsAlpha,
                            Rcpp::Named("GS.beta") = gsBeta,
                            Rcpp::Named("GS.eta") = gseta,
                            Rcpp::Named("GS.invTAUsq1") = gsInvTauSq1,
                            Rcpp::Named("GS.invTAUsq2") = gsInvTauSq2,
                            Rcpp::Named("GS.lambda.sq1") = gsLambdaStar1,
                            Rcpp::Named("GS.lambda.sq2") = gsLambdaStar2,
                            Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
                            Rcpp::Named("GS.SS1") = gsSS1,
                            Rcpp::Named("GS.SS2") = gsSS2,
                            Rcpp::Named("GS.Pi.Beta") = gsPiBeta,
                            Rcpp::Named("GS.Pi.Eta") = gsPiEta);
                          
}
