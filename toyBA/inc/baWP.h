#ifndef INC_BAWP_H_
#define INC_BAWP_H_

#include <vector>
#include "Eigen/Dense"

using namespace Eigen;

namespace baWP{

Vector2d computeEpos(const Vector2d & t_m, const Vector2d & t_t);

MatrixXd computeJack(const VectorXd & t_state, const int t_nP, const int t_mT);

VectorXd computeF(const VectorXd & t_meas, const VectorXd & t_state, const int t_nP, const int t_mT);

std::pair<int,double> gnOpt(std::vector<Vector2d> & t_points, std::vector<frame> & t_frames, int t_nIt = N_IT_MAX);
}

#endif /* INC_BAWP_H_ */
