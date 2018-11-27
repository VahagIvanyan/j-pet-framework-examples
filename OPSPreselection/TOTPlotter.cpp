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
 *  @file TOTPlotter.cpp
 */

#include <iostream>
#include <JPetWriter/JPetWriter.h>
#include "TOTPlotter.h"
#include <JPetOptionsTools/JPetOptionsTools.h>

using namespace jpet_options_tools;

using namespace std;

TOTPlotter::TOTPlotter(const char* name): JPetUserTask(name) {}

bool TOTPlotter::init()
{

  INFO("TOT plotting started");

  fOutputEvents = new JPetTimeWindow("JPetHit");

  for(int i=1; i <= getParamBank().getBarrelSlotsSize(); ++i){

    getStatistics().createHistogram(
				    new TH1F(Form("tot_strip_%d", getParamBank().getBarrelSlot(i).getID()), "TOT uncalibrated; TOT [ns]; counts", 1000., 0., 100.)
				    );

    getStatistics().createHistogram(
				    new TH1F(Form("edep_strip_%d", getParamBank().getBarrelSlot(i).getID()), "TOT uncalibrated; TOT [ns]; counts", 1000., 0., 1500.)
				    );
  }

  return true;
}

bool TOTPlotter::exec()
{
  auto timeWindow = dynamic_cast<const JPetTimeWindow*>(fEvent);

  uint n = timeWindow->getNumberOfEvents();

  for(uint i=0;i<n;++i){

    const JPetHit & hit =  dynamic_cast<const JPetHit&>(timeWindow->operator[](i));

    double tot = (hit.getSignalA().getPhe() + hit.getSignalB().getPhe());
    double edep = exp((tot*1000.+1.1483e5) / 23144.);

    getStatistics().getHisto1D(Form("tot_strip_%d", hit.getBarrelSlot().getID()))->Fill(tot);
    getStatistics().getHisto1D(Form("edep_strip_%d", hit.getBarrelSlot().getID()))->Fill(edep);

    fOutputEvents->add<JPetHit>(hit);
  }


  return true;
}


bool TOTPlotter::terminate()
{
  INFO("Done plotting TOT-s");
  return true;
}

