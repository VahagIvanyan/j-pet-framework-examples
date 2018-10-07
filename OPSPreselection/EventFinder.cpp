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
 *  @file EventFinder.cpp
 */

#include <iostream>
#include <JPetWriter/JPetWriter.h>
#include "EventFinder.h"
#include <JPetOptionsTools/JPetOptionsTools.h>

using namespace jpet_options_tools;

using namespace std;

EventFinder::EventFinder(const char* name): JPetUserTask(name) {}

bool EventFinder::init()
{

  INFO("Event finding started.");

  fOutputEvents = new JPetTimeWindow("JPetEvent");

  if (isOptionSet(fParams.getOptions(), fEventTimeParamKey)){
    kEventTimeWindow = getOptionAsFloat(fParams.getOptions(), fEventTimeParamKey);
  }
  
  getStatistics().createHistogram(
				  new TH1F("hits_per_event", "Number of Hits in Event", 20, 0.5, 20.5)
				  );

  getStatistics().createHistogram(
				  new TH1F("hits_per_event","Number of Hits in Event",20, 0.5, 20.5)
				  );

  getStatistics().createHistogram(new TH1F("tot4_3+hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_2hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_1hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_any_hits", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
	
  getStatistics().createHistogram(new TH2F("threshold_fraction_3+hits", "Fraction of fired thresholds;side A;side B", 5, -0.125, 1.125, 5, -0.125, 1.125));
  getStatistics().createHistogram(new TH2F("threshold_fraction_2hits", "Fraction of fired thresholds;side A;side B", 5, -0.125, 1.125, 5, -0.125, 1.125));
  getStatistics().createHistogram(new TH2F("threshold_fraction_1hits", "Fraction of fired thresholds;side A;side B", 5, -0.125, 1.125, 5, -0.125, 1.125));
  getStatistics().createHistogram(new TH2F("threshold_fraction_any_hits", "Fraction of fired thresholds;side A;side B", 5, -0.125, 1.125, 5, -0.125, 1.125));
    
    
    return true;
}

bool EventFinder::exec()
{

  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    //uint n = timeWindow->getNumberOfEvents();
    // for(uint i=0;i<n;++i){
    //   fHitVector.push_back(dynamic_cast<const JPetHit&>(timeWindow->operator[](i)));
    // }

    vector<JPetEvent> events = buildEvents(*timeWindow);

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


bool EventFinder::terminate()
{
  INFO("Event fiding ended.");
  return true;
}

vector<JPetEvent> EventFinder::buildEvents(const JPetTimeWindow& hits)
{

  vector<JPetEvent> eventVec;
  //	sort(hitVec.begin(), hitVec.end(), sortByTimeValue);

  int s = 0;
  int nhits = hits.getNumberOfEvents();

  while ( s < nhits ) {

    JPetEvent event;
    event.setEventType(JPetEventType::kUnknown);

    const JPetHit& startHit = dynamic_cast<const JPetHit&>(hits[s]);

    event.addHit(startHit);

    int k = 1;
    while ( s + k < nhits ) {
      const JPetHit& currentHit = dynamic_cast<const JPetHit&>(hits[s + k]);
      if (fabs(currentHit.getTime() - startHit.getTime()) < kEventTimeWindow) {
        event.addHit(currentHit);
        k++;
      } else {
        break;
      }
    }

    s += k;

    if (fSaveControlHistos) getStatistics()
      .getHisto1D("hits_per_event")
      ->Fill(event.getHits().size());

    eventVec.push_back(event);
  }

  return eventVec;
}

void EventFinder::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto & event : events) {

    // Filling of histograms
    for(auto & hit: event.getHits()){
      double tot = (hit.getSignalA().getPhe() + hit.getSignalB().getPhe());

      // fill histos for all events
      getStatistics().getHisto1D("tot4_any_hits")->Fill(tot);
      getStatistics().
	getHisto2D("threshold_fraction_any_hits")->Fill(hit.getSignalA().getQualityOfPhe(),
						       hit.getSignalB().getQualityOfPhe()
						       );

      // fill histos only 1-hit events
      if(event.getHits().size()==1){
	getStatistics().getHisto1D("tot4_1hits")->Fill(tot);
	getStatistics().
	  getHisto2D("threshold_fraction_1hits")->Fill(hit.getSignalA().getQualityOfPhe(),
						      hit.getSignalB().getQualityOfPhe()
						      );
      }

      // fill histos only for 2-hit events
      if(event.getHits().size()==2){
	getStatistics().getHisto1D("tot4_2hits")->Fill(tot);
	getStatistics().
	  getHisto2D("threshold_fraction_2hits")->Fill(hit.getSignalA().getQualityOfPhe(),
						      hit.getSignalB().getQualityOfPhe()
						      );
      }

      // fill histos for 3+ hits events
      if(event.getHits().size()>=3){
	getStatistics().getHisto1D("tot4_3+hits")->Fill(tot);
	getStatistics().
	  getHisto2D("threshold_fraction_3+hits")->Fill(hit.getSignalA().getQualityOfPhe(),
						      hit.getSignalB().getQualityOfPhe()
						      );
      }
    }

    // store events with 3+ hits
    if( event.getHits().size()>=3){
      fOutputEvents->add<JPetEvent>(event);
    }
    
  }

}
