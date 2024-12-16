//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#ifndef MUPIPROJECT_RAZZLED_H
#define MUPIPROJECT_RAZZLED_H


class Razzled {
public:
  Razzled()
  : photon_score(-1)
  , electron_score(-1)
  , proton_score(-1)
  , pion_score(-1)
  , muon_score(-1)
  , razzled_pdg(-1)
  {}

    double photon_score;
    double electron_score;
    double proton_score;
    double pion_score;
    double muon_score;
    int razzled_pdg;
};

#endif //MUPIPROJECT_RAZZLED_H
