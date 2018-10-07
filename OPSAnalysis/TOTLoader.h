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
 *  @file TOTLoader.h
 */

#ifndef TOTLOADER_H 
#define TOTLOADER_H

#include <vector>
#include <map>
#include <string>
#include <JPetUserTask/JPetUserTask.h>
#include <JPetHit/JPetHit.h>
#include <JPetEvent/JPetEvent.h>

class JPetWriter;

#ifdef __CINT__
#	define override
#endif

class TOTLoader : public JPetUserTask{
public:
  TOTLoader(const char * name);
  virtual ~TOTLoader(){}
  virtual bool init() override;
  virtual bool exec() override;
  virtual bool terminate() override;
  void fillTOThistos(const JPetEvent & event);
protected:
  const std::string fTOTnormFileParamKey = "TOTLoader_TOTnormFile_std::string";  
  std::string fTOTnormFile;
  std::map<int, double> fTOTnormFactors;
};
#endif /*  !TOTLOADER_H */










