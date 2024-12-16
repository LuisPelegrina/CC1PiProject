//
// Created by Luis Pelegrina Gutiérrez on 19/3/24.
//

#ifndef MUPIPROJECT_CHI2_H
#define MUPIPROJECT_CHI2_H


class Chi2 {
public:
  Chi2()
    :proton_score(-1)
    ,kaon_score(-1)
    ,pion_score(-1)
    ,muon_score(-1)
    {}

    double proton_score;
    double kaon_score;
    double pion_score;
    double muon_score;
};



#endif //MUPIPROJECT_CHI2_H
