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
 *  @file TOTLoader.cpp
 */

#include <iostream>
#include <JPetWriter/JPetWriter.h>
#include "TOTLoader.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetEvent/JPetEvent.h>
#include <cmath>
using namespace jpet_options_tools;

using namespace std;


double Edep(double tot){
  return exp( (tot*1000. + 1.1483e5) / 23144.);
}

double TOT(double Edep){
  return (-1.1483e5 + 23144.*log(Edep))/1000.;
}

TOTLoader::TOTLoader(const char* name): JPetUserTask(name) {}

bool TOTLoader::init()
{

  INFO("Starting application of TOT normalization factors.");

  fOutputEvents = new JPetTimeWindow("JPetEvent");

  if (isOptionSet(fParams.getOptions(),fTOTnormFileParamKey)) {
    fTOTnormFile = getOptionAsString(fParams.getOptions(), fTOTnormFileParamKey);
  }else{
    WARNING("No path to the file with TOT normalization factors was provided in user options.");
  }

  std::ifstream in_file;
  in_file.open(fTOTnormFile, std::ios::in);

  int strip;
  double corr;
  while(!in_file.eof()){
    in_file>>strip>>corr;
    fTOTnormFactors[strip] = corr;
  }
  
  in_file.close();

  // book histograms
  getStatistics().createHistogram(new TH1F("tot4_3+hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_2hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_1hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_any_hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  
  return true;
}

bool TOTLoader::exec()
{
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    uint n = timeWindow->getNumberOfEvents();

    for(uint i=0;i<n;++i){
      
      const JPetEvent & event = timeWindow->getEvent<JPetEvent>(i);

      JPetEvent new_event = event;

      std::vector<JPetHit> hits = event.getHits(); 
      
      for(auto & hit: hits){

	double tot = (hit.getSignalA().getPhe() + hit.getSignalB().getPhe());
	double edep = Edep(tot);

	edep *= fTOTnormFactors[hit.getBarrelSlot().getID()];
	double new_tot = TOT(edep);

	hit.setEnergy(new_tot);
      }

      new_event.setHits(hits);

      fillTOThistos(new_event);
      
      fOutputEvents->add<JPetEvent>(new_event);

    }
      
  } else {
    return false;
  }
  return true;
}

bool TOTLoader::terminate()
{
  INFO("Application of TOT normalization factors done.");
  return true;
}

void TOTLoader::fillTOThistos(const JPetEvent & event) {

    // Filling of histograms
    for(auto & hit: event.getHits()){
      double tot = hit.getEnergy();

      // fill histos for all events
      getStatistics().getHisto1D("tot4_any_hits")->Fill(tot);

      // fill histos only 1-hit events
      if(event.getHits().size()==1){
	getStatistics().getHisto1D("tot4_1hits")->Fill(tot);
      }

      // fill histos only for 2-hit events
      if(event.getHits().size()==2){
	getStatistics().getHisto1D("tot4_2hits")->Fill(tot);
      }

      // fill histos for 3+ hits events
      if(event.getHits().size()>=3){
	getStatistics().getHisto1D("tot4_3+hits")->Fill(tot);
      }
    }
}
