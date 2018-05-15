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
 *  @file SinogramCreatorTools.cpp
 */

#include "SinogramCreatorTools.h"

unsigned int SinogramCreatorTools::roundToNearesMultiplicity(float numberToRound, float accuracy)
{
  return std::floor((numberToRound / accuracy) + (accuracy / 2));
}

int SinogramCreatorTools::calcuateAngle(float firstX, float secondX, float firstY, float secondY, float distance)
{
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

  int angleRound = std::round(angle);
  return angleRound >= 180 ? angleRound - 180 : angleRound;
}