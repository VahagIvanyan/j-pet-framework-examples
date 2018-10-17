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
 *  @file SignalTransformer.cpp
 */

#include "JPetWriter/JPetWriter.h"
#include "SignalTransformer.h"

SignalTransformer::SignalTransformer(const char* name): JPetUserTask(name) {}

bool SignalTransformer::init()
{
  INFO("Signal transforming started: Raw to Reco and Phys");
  fOutputEvents = new JPetTimeWindow("JPetPhysSignal");

  // identify the threshold order
  std::map<int,int> m;
  for(auto& tc: getParamBank().getTOMBChannels()){
    m[(int)(tc.second->getThreshold())] = tc.second->getLocalChannelNumber();
  }

  for(auto& e: m){
    thr_order.push_back(e.second);
  }

  std::cout << "Order in SignalTransformer" << std::endl;
  for(int i=0;i<4;++i){
    std::cout << thr_order.at(i) << std::endl;
  }
  
  return true;
}

bool SignalTransformer::exec()
{
  if(auto & timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    uint n = timeWindow->getNumberOfEvents();
    for(uint i=0;i<n;++i){
      const JPetRawSignal & currSignal = dynamic_cast<const JPetRawSignal&>(timeWindow->operator[](i));
      // Make Reco Signal from Raw Signal
      auto recoSignal = createRecoSignal(currSignal);
      // Make Phys Signal from Reco Signal and save
      auto physSignal = createPhysSignal(recoSignal, thr_order.front());
      fOutputEvents->add<JPetPhysSignal>(physSignal);
    }
  }else return false;
  
  return true;
}

bool SignalTransformer::terminate()
{
  INFO("Signal transforming finished");
  return true;
}

/**
 * Method rewrites Raw Signal to Reco Signal. All fields set to -1.
 */
JPetRecoSignal SignalTransformer::createRecoSignal(const JPetRawSignal& rawSignal)
{
  JPetRecoSignal recoSignal;
  recoSignal.setRawSignal(rawSignal);
  recoSignal.setAmplitude(-1.0);
  recoSignal.setOffset(-1.0);
  recoSignal.setCharge(-1.0);
  recoSignal.setDelay(-1.0);
  return recoSignal;
}

/**
 * Method rewrites Reco Signal to Phys Signal.
 * Time of Signal set to time of the Signal Channel at the First Threshold.
 * This should be changed to something more resonable. Rest of the fields
 * set to -1, quality fields set to 0.
 */
JPetPhysSignal SignalTransformer::createPhysSignal(const JPetRecoSignal& recoSignal, int thrToUse)
{
  JPetPhysSignal physSignal;
  physSignal.setRecoSignal(recoSignal);
  physSignal.setPhe(-1.0);
  physSignal.setQualityOfPhe(0.0);
  auto times = recoSignal.getRawSignal().getTimesVsThresholdNumber(JPetSigCh::Leading);
  physSignal.setTime(times.at(thrToUse));
  physSignal.setQualityOfTime(0.0);

  // calculate TOT
  std::map<int,double> leadingPoints = recoSignal.getRawSignal().getTimesVsThresholdNumber(JPetSigCh::Leading);
  std::map<int,double> trailingPoints = recoSignal.getRawSignal().getTimesVsThresholdNumber(JPetSigCh::Trailing);

  std::vector<double> tots;
  int threshold_counter = 0;
  double TOTsum = 0.;

  for(int i=1;i<5;i++){
    auto leadSearch = leadingPoints.find(i);
    auto trailSearch = trailingPoints.find(i);
    if (leadSearch != leadingPoints.end()
	&& trailSearch != trailingPoints.end()){
      tots.push_back(trailSearch->second - leadSearch->second);
      TOTsum += (trailSearch->second - leadSearch->second) / 1000.; // in ns
      threshold_counter++;
    }
  }
  
  if(tots.size()!=0){
    physSignal.setPhe(TOTsum);
    physSignal.setQualityOfPhe(threshold_counter/4.0);
  }else{
    physSignal.setPhe(0.);
    physSignal.setQualityOfPhe(0.);
  }
  
  return physSignal;
}
