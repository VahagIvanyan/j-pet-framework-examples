#include "reconstructor.h"
#include <iostream>
#include <cassert>
#include <TMinuit.h>
#include <TMath.h>

Reconstructor * Reconstructor::_this = NULL;

Reconstructor * Reconstructor::GetInstance(){ 
  if( _this == NULL ){
    _this = new Reconstructor();
  }

  return _this; 
}


// global function for the "Old" version of Minuit
void gfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  Reconstructor * rec = Reconstructor::GetInstance();
  f = rec->fcn( par );
  
}

int Reconstructor::kinematicFit(float & chisq){
  
  int status = 0;
  
  // n_fit_params = 12;
  
  // TMinuit *gMinuit = new TMinuit(n_fit_params);
  // gMinuit->SetFCN(gfcn);
  // gMinuit->SetPrintLevel(-1);
  // gMinuit->SetErrorDef( 1.0 );
  
  // double fStripWidth = 0.7;
  // double fStripHeight = 1.9;
    
  
  // // set parameters
  // for(int i=0;i<3;i++){ // 3 gammas
  //   p0[4*i] = 0.;
  //   p0[4*i+1] = 0.;
  //   p0[4*i+2] = 0.;
  //   p0[4*i+3] = 0.;
  //   sigmas[4*i] = 0.3*fStripWidth/(2*sqrt(3));
  //   sigmas[4*i+1] = 0.3*fStripHeight/(2*sqrt(3));
  //   sigmas[4*i+2] = 0.93;
  //   sigmas[4*i+3] = 0.08;

  //   double center = 0.0;
    
  //   gMinuit->mnparm(4*i, Form("x%d", i), center, sigmas[4*i], -0.5*fStripWidth, 0.5*fStripWidth, status);
  //   gMinuit->mnparm(4*i+1, Form("y%d", i), center, sigmas[4*i+1], -0.5*fStripHeight, 0.5*fStripHeight, status);
  //   /* gMinuit->mnparm(4*i+2, Form("z%d", i), center, sigmas[4*i+2], -3.*sigmas[4*i+2], 3.*sigmas[4*i+2], status); */
  //   /* gMinuit->mnparm(4*i+3, Form("t%d", i), center, sigmas[4*i+3], -3.*sigmas[4*i+3], 3.*sigmas[4*i+3], status); */
  //   gMinuit->mnparm(4*i+2, Form("z%d", i), center, sigmas[4*i+2], 0., 0., status);
  //   gMinuit->mnparm(4*i+3, Form("t%d", i), center, sigmas[4*i+3], 0., 0., status);

  // }

  // /* gMinuit->FixParameter(3); */
  // /* gMinuit->FixParameter(7); */
  // /* gMinuit->FixParameter(11); */
  
  // Double_t arglist[10];
  // arglist[0] = 2;
  // gMinuit->mnexcm( "SET STR", arglist, 2, status );

  // arglist[0] = 10000.;
  // arglist[1] = 0.1;
  // gMinuit->mnexcm("MINIMIZE", arglist, 2, status);

  // //
  // /* for(int i=0;i>n_fit_params;++i){ */
  // /*   gMinuit->FixParameter(i); */
  // /* } */
  // /* gMinuit->Release(3); */
  // /* gMinuit->Release(7); */
  // /* gMinuit->Release(11); */

  // /* gMinuit->mnexcm("MINIMIZE", arglist, 2, status); */
  // //
  
  // Double_t amin,edm,errdef;
  // Int_t nvpar,nparx,icstat;
  // gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

  // // *minval = amin;
  
  // // get resulting parameters
  // TString name;
  // Double_t val, err, llim, ulim;
  // Int_t  aaa;
  // double *x = new double[n_fit_params];
  // for(int i=0;i<n_fit_params;i++){
  //   gMinuit->mnpout( i, name, val, err, llim, ulim, aaa );
  //   x[i] = val;
  // }

  // chisq = getChisq( x );

  // setParameters( x );


  


  
  // delete[] x;
  // delete gMinuit;
  
  return status;
  
}

double Reconstructor::getChisq(const double * x) const{


  double chi = 0.0;
  for(int i=0;i<n_fit_params;++i){
    chi += pow( (x[i]-p0[i])/sigmas[i] , 2.);
  }

  return chi;
}

double Reconstructor::fcn(const double * x){

  double lambda = 100.0;
  double constraint = 0.0;

  setParameters( x );
  reconstruct();
  
  double X = sol_hit[1].X();
  double Y = sol_hit[1].Y();
  double Z = sol_hit[1].Z();
  
  for(int i=0;i<3;++i){
    constraint += pow( sqrt(X*X + Y*Y) - fChamberRadius, 2. ) ;/// 4.0;
    /* if( Z < -5.0 ){ */
    /*   constraint += pow( Z + 5.0, 2. ) ; */
    /* }else if( Z > 5.0 ){ */
    /*   constraint += pow( Z - 5.0, 2. ) ; */
    /* } */
    
  }

  return getChisq( x ) + lambda * constraint;
  
}

void Reconstructor::setParameters(const double * x){
  for(int i=0;i<3;i++){ // 3 gammas
    gammas_shift[i].SetXYZ(x[4*i+0], x[4*i+1], x[4*i+2]);
    gammas_time_shift[i] = x[4*i+3];
  }
}

/**
 * Returns the absolute position of the corrected hit
 * in the Barrel frame of reference
 */
TVector3 Reconstructor::getGammaHit(int i) const {
  
  double angle = gammas_polar_angle[i] * TMath::DegToRad();
      
  TVector3 shift_rotated = gammas_shift[i];
  shift_rotated.RotateZ( angle );

  return gammas_hit[i] + shift_rotated;
}

Double_t Reconstructor::getGammaTime(int i) const {
  return gammas_time[i] + gammas_time_shift[i];
}


void Reconstructor::fillParamsArray(Float_t * params){
  
  for(int i=0;i<3;i++){ // 3 gammas
    params[4*i] = getGammaHit(i).X();
    params[4*i+1] = getGammaHit(i).Y();
    params[4*i+2] = getGammaHit(i).Z();
    params[4*i+3] = getGammaTime(i);
  }
}

void Reconstructor::fillShiftsArray(Float_t * shifts){
  
  for(int i=0;i<3;i++){ // 3 gammas
    shifts[4*i] = gammas_shift[i].X();
    shifts[4*i+1] = gammas_shift[i].Y();
    shifts[4*i+2] = gammas_shift[i].Z();
    shifts[4*i+3] = gammas_time_shift[i];
  }
}


// end of stuff for the kinematic fit



const Double_t Reconstructor::cvel = 29.9792458; // cm/ns

/**
 * @brief Set the length of the barrel in centimeters
 *
 * The barrel length must be set for the reconstruction to work corrrectly
 */
void Reconstructor::setBarrelLength(double length){
  fStripLength = length;
}

/**
 * @brief Set the transverse dimensions [cm] of the scintillator strips
 *
 * @param width  dimension of the strip cross-section perpendicular to the barrel radius
 * @param height dimension of the strip cross-section parallel to the barrel radius
 *
 * These dimensions must be set for the kinematic fit to work correctly.
 * General o-Ps->3gamma reconstruction does not use these values.
 */
void Reconstructor::setStripDimensions(double width, double height){
  fStripWidth = width;
  fStripHeight = height;
}

void Reconstructor::setGamma(ushort i, const JPetHit & hit){
  assert(i <= 2);
  gammas_hit[i] = hit.getPos();
  gammas_time[i] = hit.getTime() / 1000.;
  gammas_polar_angle[i] = hit.getBarrelSlot().getTheta();
  // initialize structures for the kinematic fit
  gammas_time_shift[i] = 0;
  gammas_shift[i].SetXYZ(0, 0, 0);
}

void Reconstructor::setPromptGamma(const JPetHit & hit){
  prompt_gamma_hit = hit.getPos();
  prompt_gamma_time = hit.getTime();
}

int Reconstructor::getSolution(TVector3 & point, Double_t & time, int whichSolution){
  int result = reconstruct();

  point = sol_hit[whichSolution];
  time = sol_time[whichSolution];

  return result;
};

int Reconstructor::reconstruct(){
  
  int errFlag = 0;

  // find the decay plane
  TVector3 normal = ((getGammaHit(1)-getGammaHit(0)).Cross( getGammaHit(2)-getGammaHit(0) )).Unit();
  
  // prepare transformation to the decay plane
  TVector3 z(0.,0.,1.);
  TVector3 rotAxis = normal.Cross( z );
  double angle = z.Angle( normal ); // radians
  
  TRotation rot;
  rot.Rotate(angle, rotAxis);
  
  // transform gamma hits to decay plane
  TVector3 gammas2D[3];
  for(int i=0;i<3;i++){
    gammas2D[i] = rot * getGammaHit(i);
  }
  
  // solve in 2D
  int combs[][2] = {{0,1}, {1,2}, {0,2}};
  double M[3][3];
  double D[3];
  
  // fill the matrix and constants vector
  int i,j;
  for(int k=0;k<3;++k){ // k - rows
    i = combs[k][0];
    j = combs[k][1];
    M[k][0] = 2.*( gammas2D[i].X() - gammas2D[j].X() );
    M[k][1] = 2.*( gammas2D[i].Y() - gammas2D[j].Y() );
    M[k][2] = 2.*cvel*cvel*( getGammaTime(j) - getGammaTime(i) );       
    D[k] = pow(gammas2D[i].X(),2.)
      - pow(gammas2D[j].X(),2.)
      + pow(gammas2D[i].Y(),2.)
      - pow(gammas2D[j].Y(),2.)
      - cvel*cvel*pow(getGammaTime(i),2.)
      + cvel*cvel*pow(getGammaTime(j),2.);
  }

  /*
  for(int k=0;k<3;++k){ // k - rows
    std::cout << "m1 = " << M[k][0] << std::endl;
    std::cout << "m2 = " << M[k][1] << std::endl;
    std::cout << "m3 = " << M[k][2] << std::endl;
    std::cout << "D = " << D[k] << std::endl;
  }
  */
  
  // use analytical solutions: x = Ex*t+Fx, y=Ey*t+Fy
  double Ex,Ey,Fx,Fy;
  Ex = ( M[0][2]*M[1][1]-M[0][1]*M[1][2] )/( M[0][1]*M[1][0]-M[0][0]*M[1][1] );
  Fx = ( M[0][1]*D[1]-M[1][1]*D[0] )/( M[0][1]*M[1][0]-M[0][0]*M[1][1] );
  
  Ey = ( M[0][0]*M[1][2] - M[0][2]*M[1][0] )/( M[0][1]*M[1][0]-M[0][0]*M[1][1] );
  Fy = ( M[1][0]*D[0] - M[0][0]*D[1] )/( M[0][1]*M[1][0]-M[0][0]*M[1][1] );       
  
  // find t - using ready analytical solutions
  double a,b,cc,delta;
  
  a = Ex*Ex + Ey*Ey - cvel*cvel;
  b = 2.*( Ex*(Fx-gammas2D[0].X()) + Ey*(Fy-gammas2D[0].Y()) + cvel*cvel*getGammaTime(0) );
  cc = pow(Fx-gammas2D[0].X(), 2.) + pow(Fy-gammas2D[0].Y(), 2.) - cvel*cvel*pow(getGammaTime(0), 2.);
  delta = b*b - 4.*a*cc;
    if( delta < 0. ){
      errFlag = 1;
      return errFlag;
    }
    sol_time[0] = (-1.*b - sqrt(delta))/(2.*a);
    sol_time[1] = (-1.*b + sqrt(delta))/(2.*a);
    
    for(int i = 0; i<2;++i){
      TVector3 sol2Dv( Ex*sol_time[i]+Fx, Ey*sol_time[i]+Fy, gammas2D[0].Z() );
      
      // transform the solution back to 3D
      sol_hit[i] =  rot.Inverse() * sol2Dv;
      
    }
    
    // check solution 2 for reasonability
    if( errFlag == 0 ){
      if( sol_time[1] < -20000. || sol_time[1] > 20000.  ){
	errFlag = 2;
      }else if( sol_hit[1].Perp() > getGammaHit(1).Perp() ||
		fabs( sol_hit[1].Z() ) > fStripLength/2. ){
	errFlag = 3;
      }else if( TMath::IsNaN( sol_time[1] ) ||
		TMath::IsNaN( sol_hit[1].X() ) ||
		TMath::IsNaN( sol_hit[1].Y() ) ||
		TMath::IsNaN( sol_hit[1].Z() )
		){
	errFlag = 4;
      }
    }
    
  
  return errFlag;
}


double Reconstructor::getPScreationTime(){
  double prompt_path = (prompt_gamma_hit - sol_hit[1]).Mag();
  double prompt_emission_time = prompt_gamma_time - prompt_path / cvel;
  return prompt_emission_time;
}

double Reconstructor::getPSlifetime(){
  return sol_time[1] - getPScreationTime();
}

