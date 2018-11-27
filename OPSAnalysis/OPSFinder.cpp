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
 *  @file OPSFinder.cpp
 */

#include <iostream>
#include <JPetWriter/JPetWriter.h>
#include "OPSFinder.h"
#include <JPetOptionsTools/JPetOptionsTools.h>

using namespace jpet_options_tools;

using namespace std;

OPSFinder::OPSFinder(const char* name): JPetUserTask(name) {}

bool OPSFinder::init()
{
  INFO("Event finding started.");

  fOutputEvents = new JPetTimeWindow("JPetEvent");

  if (isOptionSet(fParams.getOptions(), fEventTimeParamKey)){
    kEventTimeWindow = getOptionAsFloat(fParams.getOptions(), fEventTimeParamKey);
  }else{
    ERROR("Event time window width not provided by the user!");
    return false;
  }

  for(auto key : fTOTcutKeys){
    if (isOptionSet(fParams.getOptions(), key)){
      fTOTcuts.push_back(getOptionAsFloat(fParams.getOptions(), key));
    }else{
      ERROR(Form("TOT cut value (%s) not provided by the user!", key.c_str()));
      return false;
    }
  }
  INFO(Form("Loaded TOT cut values: (%lf, %lf) and (%lf, %lf).",
	    fTOTcuts[0], fTOTcuts[1], fTOTcuts[2], fTOTcuts[3]));  

  if (isOptionSet(fParams.getOptions(), fAngleSumCutKey)){
    fAngleSumCut = getOptionAsFloat(fParams.getOptions(), fAngleSumCutKey);
  }else{
    ERROR("Angles sum cut value not provided by the user!");
    return false;
  }

  // initialize counters
  for(int i=0;i<10;++i){
    fEventCouters[i] = 0;
  }
  
  
  getStatistics().createHistogram(
				  new TH1F("hits_per_event", "Number of Hits in event", 20, 0.5, 20.5)
				  );

  getStatistics().createHistogram(
				  new TH1F("anh_hits_per_event",
					   "Number of annihilation candidate hits in Event", 10, 0.5, 10.5)
				  );

  getStatistics().createHistogram(
				  new TH1F("dex_hits_per_event",
					   "Number of deexcitation candidate hits in Event", 10, 0.5, 10.5)
				  );
  getStatistics().createHistogram(
				  new TH1F("dex_hits_per_3anh_event",
					   "Number of deexcitation candidate hits in an Event"
					   " with 3 annihilation candidates"
					   , 10, 0.5, 10.5)
				  );  

  getStatistics().createHistogram(
				  new TH1F("dex_hits_per_final_candidate",
					   "Number of deexcitation candidate hits in an Event"
					   " after all the cuts"
					   , 10, 0.5, 10.5)
				  );  

  getStatistics().createHistogram(
				  new TH2F("dex_vs_anh_hits_per_event",
					   "Number of deexcitation vs annihilation candidate"
					   " hits in an event; annihilation; deexcitation",
					   10, 0.5, 10.5, 10, 0.5, 10.5)
				  );
  

  getStatistics().createHistogram(new TH1F("tot4_3+hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_2hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_1hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_any_hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));

  getStatistics().createHistogram(new TH1F("tot4_3+hits_nocut", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_2hits_nocut", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_1hits_nocut", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_any_hits_nocut", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  
  /************************************************************************/
  /* Relative angles histogram for 3-hit events                           */
  /************************************************************************/
  getStatistics().createHistogram(
				  new TH2F("3_hit_angles",
					   "3 Hit angles difference 1-2, 2-3",
					   360, -0.5, 359.5,
					   360, -0.5, 359.5)
				  );
  getStatistics().getHisto2D("3_hit_angles")
    ->GetXaxis()->SetTitle("Smallest angle + Second smallest angle [deg]");
  getStatistics().getHisto2D("3_hit_angles")
    ->GetYaxis()->SetTitle("Second smallest angle - Smallest angle [deg]");

  // same as above, but only for passed events
  getStatistics().createHistogram(
				  new TH2F("3_hit_angles_passed",
					   "3 Hit angles difference 1-2, 2-3: passed events",
					   360, -0.5, 359.5,
					   360, -0.5, 359.5)
				  );
  getStatistics().getHisto2D("3_hit_angles_passed")
    ->GetXaxis()->SetTitle("Smallest angle + Second smallest angle [deg]");
  getStatistics().getHisto2D("3_hit_angles_passed")
    ->GetYaxis()->SetTitle("Second smallest angle - Smallest angle [deg]");  

  /************************************************************************/
  /* Theta angle diff between hits for scattering veto                    */
  /************************************************************************/
  getStatistics().createHistogram(
				  new TH1F("theta_diffs",
					   "#Delta #theta for hits in event",
					   181, -0.5, 180.5)
				  );

  getStatistics().createHistogram(
				  new TH1F("dvt",
					   "t-d/v;t-d/v [cm]",
					   200, -5., 5.)
				  );

  /************************************************************************/
  /* Studies of single-threshold hits                                     */
  /************************************************************************/
  for(int i=0;i<6;++i){
    getStatistics().createHistogram(
				    new TH1F(Form("1thr_hits_in_event_%d", i),
					     "Number of hits having 1-threshold signals in an event",
					     6, -0.5, 5.5
					     )
				    );
  }

  getStatistics().createHistogram(
				  new TH1F("1thr_number_A",
					   "Single threshold A",
					   5, -0.5, 4.5)
				  );
  getStatistics().createHistogram(
				  new TH1F("1thr_number_B",
					   "Single threshold B",
					   5, -0.5, 4.5)
				  );

  getStatistics().createHistogram(new TH2F("tot_rel_single_strip",
					   "Relative TOTs of 2 hits found in the same scintillator;"
					   "TOT 1 [ns]; TOT 2 [ns]",
					   1000, 0., 100.,
					   1000, 0., 100.)
				  );
  
  getStatistics().createHistogram(new TH2F("tot_rel_close_strips",
					   "Relative TOTs of 2 hits found in neighbouring strips;"
					   "TOT 1 [ns]; TOT 2 [ns]",
					   1000, 0., 100.,
					   1000, 0., 100.)
				  );
  
  return true;
}

bool OPSFinder::exec()
{

  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    plotTOTs(*timeWindow);
    
    std::vector<JPetEvent> events = vetoScatterings( refineEvents(*timeWindow) );

    saveEvents(events);

    fHitVector.clear();
  } else {
    return false;
  }
  return true;
}

//sorting method
bool sortByTimeValue(JPetHit hit1, JPetHit hit2)
{
  return (hit1.getTime() < hit2.getTime());
}

void OPSFinder::study1thrHits(const JPetEvent& event, int step_id){

  int count = 0;
  for(auto & hit: event.getHits()){
    if( hit.getQualityOfEnergy() > 0. && hit.getQualityOfEnergy() < 0.15 ){
      count++;
    }
  }
  getStatistics().getHisto1D(Form("1thr_hits_in_event_%d", step_id))->Fill(count);
}


bool OPSFinder::terminate()
{
  INFO("Event fiding ended.");

  INFO(Form("[Counters] Refined events:              %d", fEventCouters[0]));
  INFO(Form("[Counters] Events after dTheta>8:       %d", fEventCouters[1]));
  INFO(Form("[Counters] Events after t-d/v:          %d", fEventCouters[2]));
  INFO(Form("[Counters] Events after same scin hits: %d", fEventCouters[3]));
  INFO(Form("[Counters] Events before angle cut:     %d", fEventCouters[4]));
  INFO(Form("[Counters] Events after angle cut:      %d", fEventCouters[5]));	
  
  return true;
}

bool OPSFinder::validateHitTOT(JPetHit & hit){

  double tot = hit.getEnergy();
    
  bool annih = false;
  bool deexc = false;

  if( tot > fTOTcuts[0] && tot < fTOTcuts[1] ){
    annih = true;
  }

  if( tot > fTOTcuts[2] && tot < fTOTcuts[3] ){
    deexc = true;
  }

  if( !annih && !deexc ){
    return false;
  }

  hit.setEnergy(tot);
  hit.setQualityOfEnergy(-1.0); // convention: Qtot < 0 for bad candidate
  
  if( annih ){ // convention: Qtot < 0.5 for annihilation
    hit.setQualityOfEnergy(0.3);
  }
  
  if( deexc ){ // convention: Qtot > 0.5 for prompt
    hit.setQualityOfEnergy(0.7);
  }

  // mark hits with signals having only 1 fired threshold
  if( hit.getSignalA().getQualityOfPhe() < 0.5 || hit.getSignalB().getQualityOfPhe() < 0.5 ){

    hit.setQualityOfEnergy(0.1);

    if(hit.getSignalA().getQualityOfPhe() < 0.5){
      std::map<int, double> points = hit.getSignalA().getRecoSignal().getRawSignal().getTimesVsThresholdNumber(JPetSigCh::Leading);
      if(points.size() == 1){
	for(auto & pair: points){
	  getStatistics().getHisto1D("1thr_number_A")->Fill(pair.first);
	}
      }
    }
    if(hit.getSignalB().getQualityOfPhe() < 0.5){
      std::map<int, double> points = hit.getSignalB().getRecoSignal().getRawSignal().getTimesVsThresholdNumber(JPetSigCh::Leading);
      if(points.size() == 1){
	for(auto & pair: points){
	  getStatistics().getHisto1D("1thr_number_B")->Fill(pair.first);
	}
      }
    }

  }
  
  // additionally, reject hits with signals having only 1 fired threshold
  // if( hit.getSignalA().getQualityOfPhe() < 0.5 || hit.getSignalB().getQualityOfPhe() < 0.5 ){
  //   return false;
  // }
  
  return true;
}

bool OPSFinder::plotTOTs(const JPetTimeWindow& preEvents)
{
  int nevents = preEvents.getNumberOfEvents();
  
  for(int entry=0; entry<nevents; ++entry){
    const JPetEvent& event = dynamic_cast<const JPetEvent&>(preEvents[entry]);
    
    const auto & hits = event.getHits();
    
    int s = 0;
    int nhits = event.getHits().size();
    
    while ( s < nhits ) {

      JPetEvent newEvent;
      newEvent.setEventType(JPetEventType::kUnknown);
      
      JPetHit startHit = hits[s];
      newEvent.addHit(startHit);
      
      int k = 1;
      while ( s + k < nhits ) {
	JPetHit currentHit = hits[s + k];
      
	if (fabs(currentHit.getTime() - startHit.getTime()) < kEventTimeWindow) {
	  newEvent.addHit(currentHit);
	  k++;
	} else {
	  break;
	}
      }

      s += k;

      fillTOThistos(newEvent);
    }
  }

  return true;
}


std::vector<JPetEvent> OPSFinder::refineEvents(const JPetTimeWindow& preEvents)
{

  vector<JPetEvent> newEventVec;

  int nevents = preEvents.getNumberOfEvents();
  
  for(int entry=0; entry<nevents; ++entry){
    const JPetEvent& event = dynamic_cast<const JPetEvent&>(preEvents[entry]);
    
    const auto & hits = event.getHits();
    
    int s = 0;
    int nhits = event.getHits().size();
    
    while ( s < nhits ) {

      JPetEvent newEvent;
      newEvent.setEventType(JPetEventType::kUnknown);

      JPetHit startHit = hits[s];

      if( validateHitTOT(startHit) ){
	newEvent.addHit(startHit);
      }else{
	s++;
	continue;
      }
      
    
      int k = 1;
      while ( s + k < nhits ) {
	JPetHit currentHit = hits[s + k];
      
	if( !validateHitTOT(currentHit) ){ // skip hits which are not annih nor prompt candidates
	  k++;
	  continue;
	}      

	if (fabs(currentHit.getTime() - startHit.getTime()) < kEventTimeWindow) {
	  newEvent.addHit(currentHit);
	  k++;
	} else {
	  break;
	}
      }

      s += k;

      fEventCouters[0]++;
      study1thrHits(newEvent, 0);
      
      newEventVec.push_back(newEvent);
    }
  }

  return newEventVec;
}

  void OPSFinder::saveEvents(const std::vector<JPetEvent>& events)
{
  int n_annih = 0;
  int n_prompt = 0;

  for (const auto & event : events) {

    fEventCouters[3]++;
    study1thrHits(event, 3);
    
    n_annih = 0;
    n_prompt = 0;
    
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

      // count annihilation and prompt candidate hits in a single event
      float hit_class = hit.getQualityOfEnergy();
      if( hit_class > 0 && hit_class < 0.5 ){
	n_annih++;
      }
      
      if( hit_class > 0.5 && hit_class < 1.0 ){
	n_prompt++;
      }

    }

    // fil number of hits of particular types
    getStatistics().getHisto1D("anh_hits_per_event")->Fill(n_annih);
    getStatistics().getHisto1D("dex_hits_per_event")->Fill(n_prompt);
    getStatistics().getHisto2D("dex_vs_anh_hits_per_event")->Fill(n_annih, n_prompt);

    if( n_annih == 3 ) {
      getStatistics().getHisto1D("dex_hits_per_3anh_event")->Fill(n_prompt);
    }
    
    // store events with exactly 3 annihilation candidates
    // and no more than one deexcitation candidate
    if( n_annih == 3 && n_prompt <= 1){

      fEventCouters[4]++;
      study1thrHits(event, 4);
      
      // study angles for 3-hit events
      //      if( analyseThreeHitEvent(event) ){

	fEventCouters[5]++;
	study1thrHits(event, 5);
	
	getStatistics().getHisto1D("dex_hits_per_final_candidate")->Fill(n_prompt);	

	fOutputEvents->add<JPetEvent>(event);
	//      }
    }
    
  }

}


bool OPSFinder::analyseThreeHitEvent(const JPetEvent& event){

  vector<JPetEvent> eventVector2;
  vector<JPetHit> hits = event.getHits();

  if(hits.size() > 3){
    
    hits.erase(
	       remove_if(hits.begin(), hits.end(),
			 [](const JPetHit& hit){
			   if( hit.getQualityOfEnergy() > 0.5 ){
			     return true;
			   }
			   return false;
			 }),
	       hits.end());
    
  } 
  
  JPetHit firstHit = event.getHits().at(0);
  JPetHit secondHit = event.getHits().at(1);
  JPetHit thirdHit = event.getHits().at(2);

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
  getStatistics().getHisto2D("3_hit_angles")->Fill(angles[0]+angles[1], angles[1] - angles[0]);

  // only take o-Ps -> 3g decay candidates
  if(angles[0]+angles[1] > fAngleSumCut){
    getStatistics().getHisto2D("3_hit_angles_passed")->Fill(angles[0]+angles[1], angles[1] - angles[0]);
    return true;
  }

  return false;
}


std::vector<JPetEvent> OPSFinder::vetoScatterings(const std::vector<JPetEvent>& events){

  std::vector<JPetEvent> returnEvents;
  
  for(auto & event : events){

    auto hits = event.getHits();
    std::set<std::vector<JPetHit>::iterator> to_erase;
    
    int n = hits.size();

    bool skip_event = false;

    double min_dvt = -10000.;
    
    for(auto i=hits.begin(); i!=hits.end();i++){
      for(auto j=i+1; j!=hits.end();j++){
    
	auto & hit1 = *i;
	auto & hit2 = *j;

	uint layer_difference = abs(hit1.getBarrelSlot().getLayer().getID() -
				    hit2.getBarrelSlot().getLayer().getID());
	
	// reject hits in the same scintillators
	if( layer_difference == 0 && hit1.getBarrelSlot() == hit2.getBarrelSlot() ){

	  // study the relative TOT of such hits on the same module
	  getStatistics().getHisto2D("tot_rel_single_strip")->Fill(hit1.getEnergy(), hit2.getEnergy());
	  
	  to_erase.insert(i);
	  to_erase.insert(j);
	  break;
	  
	}

	// reject annihilation candidate hits with close angles
	if( hit1.getQualityOfEnergy() < 0.5 && hit2.getQualityOfEnergy() < 0.5 ){
	  double d_theta = fabs(hit1.getBarrelSlot().getTheta() - hit2.getBarrelSlot().getTheta());
	  if( d_theta > 180 ){
	    d_theta = 360.0 - d_theta;
	  }
	  getStatistics().getHisto1D("theta_diffs")->Fill(d_theta);

	  if( d_theta < 8.0 ){

	    // study the relative TOT of such close hits
	  getStatistics().getHisto2D("tot_rel_close_strips")->Fill(hit1.getEnergy(), hit2.getEnergy());
	    
	    skip_event = true;
	  }
	  
	  // also try to identify scatterings based on the d - vt 
	  double d = (hit1.getPos() - hit2.getPos()).Mag();
	  double dt = fabs(hit1.getTime() - hit2.getTime()) / 1000.;
	  double dvt = dt - d/kSpeedOfLight;

	  if( fabs(dvt) < fabs(min_dvt) ){
	    min_dvt = dvt;
	  }
	  	 
	  /*
	  if( d_theta < 8.0  || dvt < 30.0 ){
	    // erase the later hit
	    if( hit1.getTime() < hit2.getTime() ){
	      to_erase.insert(j);
	    }else{
	      to_erase.insert(i);
	    }
	  }
	  */
	}
	
      } 
    } // end quadratic loop over hits

    if( !skip_event ){
      getStatistics().getHisto1D("dvt")->Fill(min_dvt);	    

      fEventCouters[1]++;
      study1thrHits(event, 1);
      
      if( min_dvt > -1.8 ){
	skip_event = true;
      }

    }
    
    if( !skip_event ){

      fEventCouters[2]++;
      
      // actually remove the vetoed hits
      for(const auto & it: to_erase){
        hits.erase(it);
      }

      JPetEvent new_event = event;
      new_event.setHits(hits);
      returnEvents.push_back(new_event);
    }
      
  } // end loop over events
  
  return returnEvents;
}


void OPSFinder::fillTOThistos(const JPetEvent & event) {

    // Filling of histograms
    for(auto & hit: event.getHits()){
      double tot = hit.getEnergy();
      
      // fill histos for all events
      getStatistics().getHisto1D("tot4_any_hits_nocut")->Fill(tot);

      // fill histos only 1-hit events
      if(event.getHits().size()==1){
	getStatistics().getHisto1D("tot4_1hits_nocut")->Fill(tot);
      }

      // fill histos only for 2-hit events
      if(event.getHits().size()==2){
	getStatistics().getHisto1D("tot4_2hits_nocut")->Fill(tot);
      }

      // fill histos for 3+ hits events
      if(event.getHits().size()>=3){
	getStatistics().getHisto1D("tot4_3+hits_nocut")->Fill(tot);
      }
    }
}
