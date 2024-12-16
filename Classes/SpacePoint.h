//
// Created by Luis Pelegrina Gutiérrez on 17/9/24.
//

#ifndef CC1PIPROJECT_SPACEPOINT_H
#define CC1PIPROJECT_SPACEPOINT_H


class SpacePoint {

  public:
    SpacePoint()
      :x(-1)
      ,y(-1)
      ,z(-1)
      ,integral(-1)
      ,sigma_integral(-1)
      ,associated_pfp_ID(-1)
      {}

    double x;
    double y;
    double z;
    double integral;
    double sigma_integral;
    int associated_pfp_ID;

    void set_space_point(double x, double y, double z, double integral, double sigma_integral, double associated_pfp_ID);
};


#endif //CC1PIPROJECT_SPACEPOINT_H
