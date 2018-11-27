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
 *  @file OPSReconstructor.h
 */

#ifndef OPSANALYZER_H 
#define OPSANALYZER_H

#include <vector>
#include <map>
#include <JPetUserTask/JPetUserTask.h>
#include <JPetHit/JPetHit.h>
#include <JPetEvent/JPetEvent.h>
#include "reconstructor.h"

class JPetWriter;

#ifdef __CINT__
#	define override
#endif

class OPSAnalyzer : public JPetUserTask{
public:
  OPSAnalyzer(const char * name);
  virtual ~OPSAnalyzer(){}
  virtual bool init() override;
  virtual bool exec() override;
  virtual bool terminate() override;
  double calcThetaSum(const std::vector<JPetHit>& hits);
 protected:
  
};
#endif /*  !OPSANALYZER_H */










