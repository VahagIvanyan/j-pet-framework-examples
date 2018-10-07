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
#include <JPetWriter/JPetWriter.h>
#include "OPSAnalyzer.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetOpsEvent/JPetOpsEvent.h>
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

    } // end loop over o-Ps events in a time window
      
  } else {
    return false;
  }
  return true;
}

bool OPSAnalyzer::terminate()
{
  INFO("Analysis of o-Ps->3g decays done.");
  return true;
}

