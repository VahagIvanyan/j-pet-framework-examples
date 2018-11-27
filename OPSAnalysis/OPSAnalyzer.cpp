/**
 *  @copyright Copyright 2016 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file OPSAnalyzer.cpp
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <JPetWriter/JPetWriter.h>
#include "OPSAnalyzer.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetOpsEvent/JPetOpsEvent.h>
#include <TMath.h>
using namespace jpet_options_tools;

using namespace std;

OPSAnalyzer::OPSAnalyzer(const char* name): JPetUserTask(name) {}

bool OPSAnalyzer::init()
{

  INFO("Analysis of previously reocnstructed o-Ps->3g decays started.");

  fOutputEvents = new JPetTimeWindow("JPetOpsEvent");

  // create histograms for annihilation position
  getStatistics().createHistogram(new TH2F("decay point XY",
					   "transverse position of the o-Ps->3g decay point;"
					   "X [cm]; Y [cm]",
					   100, -50., 50.,
					   100, -50., 50.
					   )
				  );

  getStatistics().createHistogram(new TH2F("decay point XZ",
					   "position of the o-Ps->3g decay point in XZ;"
					   "Z [cm]; X [cm]",
					   100, -50., 50.,
					   100, -50., 50.
					   )
				  );  

  getStatistics().createHistogram(new TH1F("nhits",
					   "number of hits in a 3g event",
					   5, -0.5, 4.5
					   )
				  );
  
  getStatistics().createHistogram(
				  new TH1F("t_dex_anh",
					   "Time between deexcitation and annihilation"
					   ";#Delta t [ns]",
					   400, -20.05, 19.95)
				  );  

  // create histogram for true event angles
  getStatistics().createHistogram(
				  new TH2F("3_hit_angles",
					   "3 Hit true angles difference",
					   360, -0.5, 359.5,
					   360, -0.5, 359.5)
				  );
  getStatistics().getHisto2D("3_hit_angles")
    ->GetXaxis()->SetTitle("Smallest angle + Second smallest angle [deg]");
  getStatistics().getHisto2D("3_hit_angles")
    ->GetYaxis()->SetTitle("Second smallest angle - Smallest angle [deg]");

  getStatistics().createHistogram(
				  new TH2F("3_hit_angles_2g_band",
					   "3 Hit true angles difference for evts from 2g band in the theta plot",
					   360, -0.5, 359.5,
					   360, -0.5, 359.5)
				  );
  getStatistics().getHisto2D("3_hit_angles_2g_band")
    ->GetXaxis()->SetTitle("Smallest angle + Second smallest angle [deg]");
  getStatistics().getHisto2D("3_hit_angles_2g_band")
    ->GetYaxis()->SetTitle("Second smallest angle - Smallest angle [deg]");

    getStatistics().createHistogram(
				  new TH2F("3_hit_angles_right_part",
					   "3 Hit angles difference for evts right of the 2g  band in the theta plot",
					   360, -0.5, 359.5,
					   360, -0.5, 359.5)
				  );
  getStatistics().getHisto2D("3_hit_angles_right_part")
    ->GetXaxis()->SetTitle("Smallest angle + Second smallest angle [deg]");
  getStatistics().getHisto2D("3_hit_angles_right_part")
    ->GetYaxis()->SetTitle("Second smallest angle - Smallest angle [deg]");
  
  return true;
}

bool OPSAnalyzer::exec()
{

  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    uint n = timeWindow->getNumberOfEvents();

    for(uint i=0;i<n;++i){
      const JPetOpsEvent & event = timeWindow->getEvent<JPetOpsEvent>(i);

      getStatistics().getHisto1D("nhits")->Fill(event.getHits().size());
      auto hits = event.getHits();

      TVector3 anh_point = event.getAnnihilationPoint();
      
      getStatistics().getHisto2D("decay point XY")->Fill(anh_point.X(),
							 anh_point.Y());

      getStatistics().getHisto2D("decay point XZ")->Fill(anh_point.Z(),
							 anh_point.X());

      // calculate event time
      if( hits.size() == 4 ){ // only for events with de-excitation
	double dt = -1000.;
	for(auto & hit: event.getHits()){
	  float hit_class = hit.getQualityOfEnergy();
	  if( hit_class > 0.5 && hit_class < 1.0 ){ // prompt gamma
	    double r = (anh_point - hit.getPos()).Mag();
	    const double c = 29.9792458; // cm  / ns
	    double t_prompt = hit.getTime() - r/c;
	    dt = event.getAnnihilationTime() - t_prompt;
	  }
      
	}
      
	getStatistics().getHisto1D("t_dex_anh")->Fill(dt);
      }

      //
      // calculate true event angles
      //
      
      // remove deexcitation photon candidate
      hits.erase(
		 remove_if(hits.begin(), hits.end(),
			   [](const JPetHit& hit){
			     if( hit.getQualityOfEnergy() > 0.5 ){
			       return true;
			     }
			     return false;
			   }),
		 hits.end());

      assert(hits.size() == 3);

      TVector3 k[3];
      for(int i=0;i<3;++i){
	k[i] = hits.at(i).getPos() - event.getAnnihilationPoint();
      }
      std::vector<double> angles(3, -1.);
      angles[0] = 180.*k[0].Angle(k[1]) / TMath::Pi();
      angles[1] = 180.*k[1].Angle(k[2]) / TMath::Pi();
      angles[2] = 180.*k[0].Angle(k[2]) / TMath::Pi();
      std::sort(angles.begin(), angles.end());
      assert(angles[0] <= angles[1] <= angles[2]);
      
      double sum_2_least = angles[0] + angles[1];
      double diff_2_least = angles[1] - angles[0];
      if(sum_2_least < 180.){
	angles[2] = 360. - angles[2];
      }
      
      getStatistics().getHisto2D("3_hit_angles")->Fill(sum_2_least, diff_2_least);

      // for comparison, calculate sum of detector thetas in the XY plane
      double theta_sum = calcThetaSum(hits);
      if( fabs(theta_sum-180.) < 34. ){
	getStatistics().getHisto2D("3_hit_angles_2g_band")->Fill(sum_2_least, diff_2_least);
      }
      if( theta_sum > 215. ){
	getStatistics().getHisto2D("3_hit_angles_right_part")->Fill(sum_2_least, diff_2_least);
      }
      
    } // end loop over o-Ps events in a time window
      
  } else {
    return false;
  }
  return true;
}

double OPSAnalyzer::calcThetaSum(const std::vector<JPetHit>& hits){

  JPetHit firstHit = hits.at(0);
  JPetHit secondHit = hits.at(1);
  JPetHit thirdHit = hits.at(2);

  std::vector<double> angles;
  angles.push_back(firstHit.getBarrelSlot().getTheta());
  angles.push_back(secondHit.getBarrelSlot().getTheta());
  angles.push_back(thirdHit.getBarrelSlot().getTheta());
  std::sort( angles.begin(), angles.begin() +3 );
  float theta_1_2 = angles[1] - angles[0];
  float theta_2_3 = angles[2] - angles[1];
  float theta_3_1 = 360 - theta_1_2 - theta_2_3;
  angles.clear();
  angles.push_back(theta_1_2);
  angles.push_back(theta_2_3);
  angles.push_back(theta_3_1);
  std::sort( angles.begin(), angles.begin() +3 );

  return angles[0]+angles[1];
}

bool OPSAnalyzer::terminate()
{
  INFO("Analysis of o-Ps->3g decays done.");
  return true;
}

