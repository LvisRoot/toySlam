/*
 * ba.h
 *
 *  Created on: Sep 5, 2018
 *      Author: l_vis
 */

#ifndef INC_BA_H_
#define INC_BA_H_


/*
 * main.cpp
 *
 *  Created on: Sep 2, 2018
 *      Author: l_vis
 */

#include <vector>

#include "Eigen/Dense"

const static double L_LM = 1e-5;

const static double E_TOL = 1e-6;
const static int N_IT_MAX = 1000;

const static double CAM_FX = 1;
const static int CAM_CX = 0;

const static double CI_COV = 0.5;



using namespace Eigen;

struct frame{
	Vector2d _t;
	std::vector<Vector2d> _tMeas;
	std::vector<double> _qMeas;
};

typedef struct frame frame;


Vector2d computePointJack(const Vector2d & t_point, const Vector2d & t_pose);

Vector2d computePoseJack(const Vector2d & t_point, const Vector2d & t_pose);

MatrixXd computeJack(const VectorXd & t_state, const int t_nP, const int t_mT);

double project(const Vector2d & t_p, const Vector2d & t_t);

double computeErep(const double & t_q, const Vector2d & t_p, const Vector2d & t_t);

VectorXd computeF(const VectorXd & t_meas, const VectorXd & t_state, const int t_nP, const int t_mT);

std::pair<int,double> gnOpt(std::vector<Vector2d> & t_points, std::vector<frame> & t_frames, int t_nIt = N_IT_MAX);




#endif /* INC_BA_H_ */
