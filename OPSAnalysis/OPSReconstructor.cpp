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
 *  @file OPSReconstructor.cpp
 */

#include <iostream>
#include <JPetWriter/JPetWriter.h>
#include "OPSReconstructor.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetOpsEvent/JPetOpsEvent.h>
using namespace jpet_options_tools;

using namespace std;

OPSReconstructor::OPSReconstructor(const char* name): JPetUserTask(name) {}

bool OPSReconstructor::init()
{

  INFO("Reconstruction of o-Ps->3g decays started.");

  fOutputEvents = new JPetTimeWindow("JPetOpsEvent");

  fReconstructor = Reconstructor::GetInstance();
  double scin_length = getParamBank().getScintillators().begin()->second->getScinSize(JPetScin::kLength);
  scin_length *= 0.1; // in cm
  fReconstructor->setBarrelLength( scin_length );
  fReconstructor->setChamberRadius( 10.0 );

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

  // XY projections for different Z
  getStatistics().createHistogram(new TH2F("xy_0",
					   "o-Ps->3g XY for Z in (-1,+1) cm;"
					   "X [cm]; Y [cm]",
					   100, -50., 50.,
					   100, -50., 50.
					   )
				  );
  
  
  for(int i=1;i<7;++i){
    
    int z0 = 1 + 2 * (i-1);
    int z1 = z0 + 2;
    getStatistics().createHistogram(new TH2F(Form("xy_+%d", i),
					     Form("o-Ps->3g XY for Z in (%d,%d) cm;"
						  "X [cm]; Y [cm]", z0, z1),
					     100, -50., 50.,
					     100, -50., 50.
					     )
				    );
    getStatistics().createHistogram(new TH2F(Form("xy_-%d", i),
					     Form("o-Ps->3g XY for Z in (-%d,-%d) cm;"
						  "X [cm]; Y [cm]", z1, z0),
					     100, -50., 50.,
					     100, -50., 50.
					     )
				    );
    
  }
  
  return true;
}

bool OPSReconstructor::exec()
{

  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    uint n = timeWindow->getNumberOfEvents();

    for(uint i=0;i<n;++i){
      const JPetEvent & event = timeWindow->getEvent<JPetEvent>(i);
      
      // sanity check
      uint n_hits = event.getHits().size();
      assert(n_hits == 3 || n_hits == 4 );

      // remove deexcitation photon candidate
      auto hits = event.getHits();
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
      
      //
      bool z_ok = true;

      std::vector<JPetHit*> annih_hits;
      
      for(auto & hit : event.getHits()){
	if( fabs(hit.getPos().Z()) > 23.0 ){
	  z_ok = false;
	}
      }

      // skip events with bad Z of hits
      if(!z_ok) continue;
      
      // setup reconstruction
      fReconstructor->setGamma(0, hits.at(0));
      fReconstructor->setGamma(1, hits.at(1));
      fReconstructor->setGamma(2, hits.at(2));

      int error = 0;
      TVector3 sol[2];
      double t[2];
      error = fReconstructor->getSolution(sol[1], t[1], 1);

      if(error == 0){
	getStatistics().getHisto2D("decay point XY")->Fill(sol[1].X(), sol[1].Y());
	getStatistics().getHisto2D("decay point XZ")->Fill(sol[1].Z(), sol[1].X());

	// fill appropriate XY projection
	int projection = (sol[1].Z()-1.0)/2.0;
	if( abs(projection) < 7 ){
	  TString histo_name;
	  if( projection == 0 ) histo_name = "xy_0";
	  else if( projection > 0) histo_name = Form("xy_+%d", projection);
	  else histo_name = Form("xy_-%d", abs(projection));
	  getStatistics().getHisto2D(histo_name)->Fill(sol[1].X(), sol[1].Y());
	}
	
	JPetOpsEvent ops_event(event);
	ops_event.setAnnihilationPoint(sol[1]);
	ops_event.setAnnihilationTime(t[1]);
	
	fOutputEvents->add<JPetOpsEvent>(ops_event);
	
      }
    } // end loop over events in a time window
      
  } else {
    return false;
  }
  return true;
}

bool OPSReconstructor::terminate()
{
  INFO("Reconstruction of o-Ps->3g decays done.");
  return true;
}

