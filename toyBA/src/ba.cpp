/*
 * main.cpp
 *
 *  Created on: Sep 2, 2018
 *      Author: l_vis
 */

#include "ba.h"

#include <vector>
#include <iostream>
#include <random>

#include <chrono>

#include "Eigen/Dense"


Vector2d computePointJack(const Vector2d & t_point, const Vector2d & t_pose){

	double dx = t_point[0] - t_pose[0];
	double dz = t_point[1] - t_pose[1];

	return Vector2d(-CAM_FX / dz, CAM_FX * dx / (dz * dz));
}

Vector2d computePoseJack(const Vector2d & t_point, const Vector2d & t_pose){

	return -computePointJack(t_point,t_pose);
}

MatrixXd computeJack(const VectorXd & t_state, const int t_nP, const int t_mT){

	MatrixXd jack(t_nP*t_mT,2*(t_nP+t_mT));

	for(int j = 0; j < t_mT ; ++j){
		for(int i = 0; i < t_nP ; ++i){

			jack.block<1,2>(t_nP*j+i,2*j) = computePoseJack(t_state.segment<2>(2*(t_mT+i)),t_state.segment<2>(2*j)).transpose();

			jack.block<1,2>(t_nP*j+i,2*(t_mT+i)) = computePointJack(t_state.segment<2>(2*(t_mT+i)),t_state.segment<2>(2*j)).transpose();
		}
	}
	return jack;
}

double project(const Vector2d & t_p, const Vector2d & t_t){
	return CAM_FX * (t_p[0] - t_t[0]) / (t_p[1] - t_t[1]) + CAM_CX;
}

double computeErep(const double & t_q, const Vector2d & t_p, const Vector2d & t_t){

	return t_q - project(t_p,t_t);
}

VectorXd computeF(const VectorXd & t_meas, const VectorXd & t_state, const int t_nP, const int t_mT){

	VectorXd eRep(t_nP*t_mT);

	for(int j = 0; j < t_mT ; ++j){
		for(int i = 0; i < t_nP ; ++i){

			eRep[t_nP*j+i] = computeErep(t_meas[t_nP*j+i],t_state.segment<2>(2*(t_mT+i)),t_state.segment<2>(2*j));
		}
	}
	return eRep;
}

std::pair<int,double> gnOpt(std::vector<Vector2d> & t_points, std::vector<frame> & t_frames, int t_nIt){

	const double 	m_xp = 0,
					m_zp = 0,
					m_xt = 0,
					m_zt = 0;
	const double 	cov_xp = CI_COV,
					cov_zp = CI_COV,
					cov_xt = CI_COV,
					cov_zt = CI_COV;
	std::default_random_engine gen;
	std::normal_distribution<double> 	xpDist(m_xp, cov_xp),
										zpDist(m_zp, cov_zp),
										xtDist(m_xt, cov_xt),
										ztDist(m_zt, cov_zt);

	int nP = t_points.size();
	int mT = t_frames.size();

	if(t_nIt > N_IT_MAX)
		t_nIt = N_IT_MAX;

	VectorXd state (2*(nP+mT));
	VectorXd meas (nP*mT);


	for (int j = 0; j < mT; ++j) {
		state.segment<2>(2*j) = t_frames[j]._t + Vector2d(xtDist(gen),ztDist(gen));
	}

	for (int i = 0; i < nP; ++i) {
		state.segment<2>(2*(mT+i)) = t_points[i] + Vector2d(xpDist(gen),zpDist(gen));
	}

	for(int j = 0; j < mT ; ++j){
		for(int i = 0; i < nP ; ++i){

			meas[nP*j+i] = t_frames[j]._qMeas[i];
		}
	}

	std::cout << "\n State CI: " << state.transpose();

	//VectorXd fErrOld(nP*mT);

	double eNorm = std::numeric_limits<double>::max();
	int nIts = 0;

	for (int i = 0; i < t_nIt && eNorm >= E_TOL; ++i) {

		VectorXd fErr = computeF(meas,state,nP,mT);
		MatrixXd jack = computeJack(state,nP,mT);
		MatrixXd hess = jack.transpose()*jack + L_LM * MatrixXd::Identity(jack.cols(),jack.cols());

		state += hess.ldlt().solve(-jack.transpose()*fErr);

		eNorm = fErr.norm();
		nIts = i+1;
	}

	for (int j = 0; j < mT; ++j)
		 t_frames[j]._t = state.segment<2>(2 * j);

	for (int i = 0; i < nP; ++i)
		t_points[i] = state.segment<2>(2 * (mT + i));

	return std::pair<int,double>(nIts,eNorm);
}
