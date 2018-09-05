#include "ba.h"
#include "baWP.h"

#include <vector>
#include <iostream>
#include <random>

#include <chrono>

#include "Eigen/Dense"


using namespace Eigen;

typedef struct frame frame;


MatrixXd baWP::computeJack(const VectorXd & t_state, const int t_nP, const int t_mT){

	MatrixXd jack(t_nP*t_mT+2,2*(t_nP+t_mT));

	for(int j = 0; j < t_mT ; ++j){
		for(int i = 0; i < t_nP ; ++i){

			jack.block<1,2>(t_nP*j+i,2*j) = computePoseJack(t_state.segment<2>(2*(t_mT+i)),t_state.segment<2>(2*j)).transpose();

			jack.block<1,2>(t_nP*j+i,2*(t_mT+i)) = computePointJack(t_state.segment<2>(2*(t_mT+i)),t_state.segment<2>(2*j)).transpose();
		}
	}

	jack.block<2,2>(t_nP*t_mT,0) = -Matrix2d::Identity();
	return jack;
}

Vector2d baWP::computeEpos(const Vector2d & t_m, const Vector2d & t_t){

	return t_m - t_t;
}

VectorXd baWP::computeF(const VectorXd & t_meas, const VectorXd & t_state, const int t_nP, const int t_mT){

	//TODO: compute size depending of n pose measurements
	VectorXd eFun(t_nP*t_mT + 2);

	for(int j = 0; j < t_mT ; ++j){
		for(int i = 0; i < t_nP ; ++i){

			eFun[t_nP*j+i] = 	computeErep(t_meas[t_nP*j+i],t_state.segment<2>(2*(t_mT+i)),t_state.segment<2>(2*j));

		}
	}

	//TODO: include multiple pose measurements
	eFun.segment<2>(t_nP*t_mT) = baWP::computeEpos(t_meas.segment<2>(t_nP*t_mT), t_state.segment<2>(0));

	return eFun;
}

std::pair<int,double> baWP::gnOpt(std::vector<Vector2d> & t_points, std::vector<frame> & t_frames, int t_nIt){

	const double 	cov_xp = CI_COV,
					cov_zp = CI_COV,
					cov_xt = CI_COV,
					cov_zt = CI_COV;
	std::default_random_engine gen;
	std::normal_distribution<double> 	xpDist(0, cov_xp),
										zpDist(0, cov_zp),
										xtDist(0, cov_xt),
										ztDist(0, cov_zt);

	int nP = t_points.size();
	int mT = t_frames.size();

	if(t_nIt > N_IT_MAX)
		t_nIt = N_IT_MAX;

	VectorXd state (2*(nP+mT));
	VectorXd meas (nP*mT+2);


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

	meas.segment<2>(nP*mT) = t_frames[0]._tMeas[0];

	std::cout << "\n State CI: " << state.transpose();

	//VectorXd fErrOld(nP*mT);

	double eNorm = std::numeric_limits<double>::max();
	int nIts = 0;

	for (int i = 0; i < t_nIt && eNorm >= E_TOL; ++i) {

		VectorXd fErr = baWP::computeF(meas,state,nP,mT);
		MatrixXd jack = baWP::computeJack(state,nP,mT);
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
