//
// Created by Luis Pelegrina Gutiérrez on 18/9/24.
//

#ifndef CC1PIPROJECT_CALOPOINT_H
#define CC1PIPROJECT_CALOPOINT_H

#include <iostream>

class CaloPoint {
  public:
    CaloPoint()
    :dEdx(-1)
    ,dQdx(-1)
    ,pitch(-1)
    ,residual_range(-1)
    ,x_pos(-1)
    ,y_pos(-1)
    ,z_pos(-1)
    {}
    double dEdx;
    double dQdx;
    double pitch;
    double residual_range;
    double x_pos;
    double y_pos;
    double z_pos;

};


#endif //CC1PIPROJECT_CALOPOINT_H
