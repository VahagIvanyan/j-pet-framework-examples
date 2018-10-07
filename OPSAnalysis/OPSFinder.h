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
 *  @file OPSFinder.h
 */

#ifndef OPSFINDER_H 
#define OPSFINDER_H 

#include <vector>
#include <map>
#include <JPetUserTask/JPetUserTask.h>
#include <JPetHit/JPetHit.h>
#include <JPetEvent/JPetEvent.h>

class JPetWriter;

#ifdef __CINT__
#	define override
#endif

class OPSFinder : public JPetUserTask{
public:
	OPSFinder(const char * name);
	virtual ~OPSFinder(){}
	virtual bool init() override;
	virtual bool exec() override;
	virtual bool terminate() override;
protected:
  double kEventTimeWindow = 5000.0; //ps
  const std::string fEventTimeParamKey = "OPSFinder_FineEventTime_float";

  const std::vector<std::string> fTOTcutKeys = {
    "OPSFinder_TOTanh_low_float",
    "OPSFinder_TOTanh_high_float",
    "OPSFinder_TOTdex_low_float",
    "OPSFinder_TOTdex_high_float"
  };

  const std::string fAngleSumCutKey = "OPSFinder_angles_sum_cut_float";
  float fAngleSumCut;
  
  std::vector<float> fTOTcuts;

  std::vector<JPetHit> fHitVector;
  bool fSaveControlHistos = true;
  void saveEvents(const std::vector<JPetEvent>& event);
  std::vector<JPetEvent> refineEvents(const JPetTimeWindow & preEvents);
  bool plotTOTs(const JPetTimeWindow & preEvents);
  bool validateHitTOT(JPetHit & hit);
  bool analyseThreeHitEvent(const JPetEvent& event);
  std::vector<JPetEvent> vetoScatterings(const std::vector<JPetEvent>& events);
  void fillTOThistos(const JPetEvent & event);
  void study1thrHits(const JPetEvent& event, int step_id);
  
  const double kSpeedOfLight = 29.9792458; // cm  / ns

  int fEventCouters[10];
  
};
#endif /*  !OPSFINDER_H */










