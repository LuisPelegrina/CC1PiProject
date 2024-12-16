//
// Created by Luis Pelegrina Gutiérrez on 10/10/24.
//

#include "RefinedSlice.h"
#ifndef CC1PIPROJECT_HOUGHALGORITHM_H
#define CC1PIPROJECT_HOUGHALGORITHM_H


class HoughAlgorithm {
public:
  vector<RegionOfInterest> create_ROI(RefinedSlice refined_slice);
};


#endif //CC1PIPROJECT_HOUGHALGORITHM_H
