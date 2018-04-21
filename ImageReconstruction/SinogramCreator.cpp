/**
 *  @copyright Copyright 2017 The J-PET Framework Authors. All rights reserved.
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
 *  @file SinogramCreator.cpp
 */

#include "SinogramCreator.h"
#include <TH2I.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1I.h>
using namespace jpet_options_tools;

SinogramCreator::SinogramCreator(const char* name) : JPetUserTask(name) {}

SinogramCreator::~SinogramCreator() {}

bool SinogramCreator::init()
{
  setUpOptions();
  fOutputEvents = new JPetTimeWindow("JPetEvent");

  getStatistics().createHistogram(new TH2I("reconstuction_histogram",
                                  "reconstuction histogram",
                                  std::ceil(fReconstructionLayerRadius * 2 * (1.f / fReconstructionDistanceAccuracy)) + 1, 0.f, fReconstructionLayerRadius,
                                  kReconstructionMaxAngle, 0, kReconstructionMaxAngle));

  getStatistics().createHistogram(new TH1F("pos_dis", "Position distance real data", (fReconstructionLayerRadius) * 10 * 5, 0.f, fReconstructionLayerRadius));
  getStatistics().createHistogram(new TH1F("angle", "Position angle real data", kReconstructionMaxAngle, 0, kReconstructionMaxAngle));

  getStatistics().getObject<TH2I>("reconstuction_histogram")->SetBit(TH2::kCanRebin);

  getStatistics().getObject<TH1F>("angle")->SetBit(TH1::kCanRebin);
  getStatistics().getObject<TH1F>("pos_dis")->SetBit(TH1::kCanRebin);

  return true;
}

bool SinogramCreator::exec()
{
  int maxDistanceNumber = std::ceil(fReconstructionLayerRadius * 2 * (1.f / fReconstructionDistanceAccuracy)) + 1;
  if (fSinogram == nullptr) {
    fSinogram = new SinogramResultType*[fZSplitNumber];
    for (int i = 0; i < fZSplitNumber; i++) {
      fSinogram[i] = new SinogramResultType(maxDistanceNumber, (std::vector<unsigned int>(kReconstructionMaxAngle, 0)));
    }
  }
  if (const auto& timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    unsigned int numberOfEventsInTimeWindow = timeWindow->getNumberOfEvents();
    for (unsigned int i = 0; i < numberOfEventsInTimeWindow; i++) {
      auto event = dynamic_cast<const JPetEvent&>(timeWindow->operator[](static_cast<int>(i)));
      auto hits = event.getHits();
      if (hits.size() == 2) {
        const auto& firstHit = hits[0];
        const auto& secondHit = hits[1];
        float firstX = firstHit.getPosX();
        float firstY = firstHit.getPosY();
        float secondX = secondHit.getPosX();
        float secondY = secondHit.getPosY(); // copy positions

        float firstZ = firstHit.getPosZ();
        float secondZ = secondHit.getPosZ();

        for (int i = 0; i < fZSplitNumber; i ++) {
          if (firstZ >= fZSplitRange[i].first && firstZ <= fZSplitRange[i].second
              && secondZ >= fZSplitRange[i].first && secondZ <= fZSplitRange[i].second ) {
            float distance = (secondX * firstY) - (secondY * firstX);

            float norm = std::sqrt(std::pow((secondY - firstY), 2) + std::pow((secondX - firstX), 2));
            if (norm > 0.f) {
              distance /= norm;
              getStatistics().getObject<TH1F>("pos_dis")->Fill(distance);
              float angle = 0.f;
              if (std::abs(secondX - firstX) > EPSILON)
                angle = std::atan((firstY - secondY) / (secondX - firstX));

              if (distance > 0.f)
                angle = angle + M_PI / 2.f;
              else
                angle = angle + 3.f * M_PI / 2.f;

              if (angle > M_PI) {
                angle = angle - M_PI;
              }
              angle *= 180.f / M_PI;
              if (angle >= 90.f)
                distance = -distance;
              getStatistics().getObject<TH1F>("angle")->Fill(angle);
              int distanceRound = SinogramCreatorTools::roundToNearesMultiplicity(distance + fReconstructionLayerRadius, fReconstructionDistanceAccuracy);
              int thetaNumber = std::round(angle);
              if (thetaNumber >= 180)
                thetaNumber -= 180;
              if (thetaNumber >= kReconstructionMaxAngle) {
                std::cout << "Angle: " << angle << " rounded angle: " << thetaNumber << " x1: " << firstX << " y1: " << firstY << " x2: " << secondX << " y2: " << secondY << std::endl;
              }
              if (distanceRound >= maxDistanceNumber) {
                std::cout << "Distance round: " << distanceRound << " distance: " << distance << " norm: " << norm << " x1: " << firstX << " y1: " << firstY << " x2: " << secondX << " y2: " << secondY << std::endl;
              }
              fCurrentValueInSinogram[i] = ++fSinogram[i]->at(distanceRound).at(thetaNumber);
              if (fCurrentValueInSinogram[i] >= fMaxValueInSinogram[i])
                fMaxValueInSinogram[i] = fCurrentValueInSinogram[i];                                                                 // save max value of sinogram
              getStatistics().getObject<TH2I>("reconstuction_histogram")->Fill(distance + fReconstructionLayerRadius, angle); //add to histogram
            }
          } else {
            continue;
          }
        }


      }
    }
  } else {
    ERROR("Returned event is not TimeWindow");
    return false;
  }
  return true;
}

bool SinogramCreator::checkLayer(const JPetHit& hit)
{
  return hit.getBarrelSlot().getLayer().getID() == 1;
}

bool SinogramCreator::terminate()
{
  for (int i = 0; i < fZSplitNumber; i ++) {
    std::ofstream res(fOutFileName + std::to_string(i) + ".ppm");
    res << "P2" << std::endl;
    res << (*fSinogram[i])[0].size() << " " << fSinogram[i]->size() << std::endl;
    res << fMaxValueInSinogram[i] << std::endl;
    for (unsigned int i = 0; i < fSinogram[i]->size(); i++) {
      for (unsigned int j = 0; j < (*fSinogram[i])[0].size(); j++) {
        res << (*fSinogram[i])[i][j] << " ";
      }
      res << std::endl;
    }
    res.close();
  }
  for (int i = 0; i < fZSplitNumber; i++) {
    delete[] fSinogram[i];
  }

  delete [] fMaxValueInSinogram;

  return true;
}

void SinogramCreator::setUpOptions()
{
  auto opts = getOptions();
  if (isOptionSet(opts, kOutFileNameKey)) {
    fOutFileName = getOptionAsString(opts, kOutFileNameKey);
  }

  if (isOptionSet(opts, kReconstructionLayerRadiusKey)) {
    fReconstructionLayerRadius = getOptionAsFloat(opts, kReconstructionLayerRadiusKey);
  }

  if (isOptionSet(opts, kReconstructionDistanceAccuracy)) {
    fReconstructionDistanceAccuracy = getOptionAsFloat(opts, kReconstructionDistanceAccuracy);
  }

  if (isOptionSet(opts, kZSplitNumber)) {
    fZSplitNumber = getOptionAsInt(opts, kZSplitNumber);
  }

  fMaxValueInSinogram = new int[fZSplitNumber];
  fCurrentValueInSinogram = new int[fZSplitNumber];
  const float maxZRange = 25.f;
  float range = 2.f * maxZRange;
  for (int i = 0; i < fZSplitNumber; i ++) {
    float rangeStart = (i * range) - maxZRange;
    float rangeEnd = ((i + 1) * range) - maxZRange;
    fZSplitRange.push_back(std::make_pair(rangeStart, rangeEnd));
    fCurrentValueInSinogram[i] = 0;
  }
}
