#include "../../Includes.h"

void distance_test()
{
  // Punto P
  TVector3 P(1.0, 0, 0);

  // Punto Q en la recta
  TVector3 Q(0.0, 0.0, 0.0);

  // Vector director de la recta
  TVector3 v(1.0, 1.0, 0.0);

  // Calcula el vector PQ
  TVector3 PQ = P - Q;

  // Calcula el producto cruzado entre PQ y v
  TVector3 crossProduct = PQ.Cross(v);

  // Calcula la distancia: norma del producto cruzado dividido por la norma del vector director
  double distance = crossProduct.Mag() / v.Mag();

  cout << distance << endl;
}