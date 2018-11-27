/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
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
 *  @file main.cpp
 */

#include <JPetManager/JPetManager.h>
#include "TOTLoader.h"
#include "OPSFinder.h"
#include "OPSReconstructor.h"
#include "OPSAnalyzer.h"

using namespace std;

int main(int argc, const char* argv[])
{
  JPetManager& manager = JPetManager::getManager();

  manager.registerTask<TOTLoader>("TOTLoader");
  manager.registerTask<OPSFinder>("OPSFinder");
  manager.registerTask<OPSReconstructor>("OPSReconstructor");
  manager.registerTask<OPSAnalyzer>("OPSAnalyzer");

  //  manager.useTask("TOTLoader", "pre.evt", "pre.evt.tot");
  manager.useTask("OPSFinder", "pre.evt", "ops.cand.evt");
  manager.useTask("OPSReconstructor", "ops.cand.evt", "ops.rec.evt");
  manager.useTask("OPSAnalyzer", "ops.rec.evt", "ops.ana.evt");
  
  manager.run(argc, argv);
}
