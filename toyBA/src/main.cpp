/*
 * main.cpp
 *
 *  Created on: Sep 2, 2018
 *      Author: l_vis
 */

#include <vector>
#include <iostream>
#include <random>

#include <chrono>

#include "Eigen/Dense"


const static double E_TOL = 1e-5;
const static int N_IT_MAX = 1000;

const static double CAM_FX = 1;
const static int CAM_CX = 0;

const static double CI_COV = 0.5;



using namespace Eigen;

struct frame{
	Vector2d _pose;
	std::vector<double> _meas;
};

typedef struct frame frame;


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

double gnOpt(std::vector<Vector2d> & t_points, std::vector<frame> & t_frames, int t_nIt = N_IT_MAX){

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

	double l_LM = 1e-5;

	if(t_nIt > N_IT_MAX)
		t_nIt = N_IT_MAX;

	VectorXd state (2*(nP+mT));
	VectorXd meas (nP*mT);


	for (int j = 0; j < mT; ++j) {
		state.segment<2>(2*j) = t_frames[j]._pose + Vector2d(xtDist(gen),ztDist(gen));
	}

	for (int i = 0; i < nP; ++i) {
		state.segment<2>(2*(mT+i)) = t_points[i] + Vector2d(xpDist(gen),zpDist(gen));
	}

	for(int j = 0; j < mT ; ++j){
		for(int i = 0; i < nP ; ++i){

			meas[nP*j+i] = t_frames[j]._meas[i];
		}
	}

	std::cout << "\n State CI: " << state.transpose();

	//VectorXd fErrOld(nP*mT);

	double eNorm = std::numeric_limits<double>::max();

	for (int i = 0; i < t_nIt && eNorm >= E_TOL; ++i) {

		VectorXd fErr = computeF(meas,state,nP,mT);
		MatrixXd jack = computeJack(state,nP,mT);
		MatrixXd hess = jack.transpose()*jack + l_LM * MatrixXd::Identity(jack.cols(),jack.cols());

		state += hess.ldlt().solve(-jack.transpose()*fErr);

		eNorm = fErr.norm();
	}

	for (int j = 0; j < mT; ++j)
		 t_frames[j]._pose = state.segment<2>(2 * j);

	for (int i = 0; i < nP; ++i)
		t_points[i] = state.segment<2>(2 * (mT + i));

	return eNorm;
}


int main(int argc, char **argv) {

	int nP = atoi(argv[1]);
	int mT = atoi(argv[2]);
	int nIts = atoi(argv[3]);

	std::vector<Vector2d> points, poses;
	std::vector<frame> frames;

    const double 	m_xp = 4,
					m_zp = 3,
    				m_xt = 3,
					m_zt = -2;
    const double 	cov_xp = 2,
					cov_zp = 0.1,
					cov_xt = 3,
					cov_zt = 0.1;
    std::default_random_engine gen;
    std::normal_distribution<double> 	xpDist(m_xp, cov_xp),
										zpDist(m_zp, cov_zp),
										xtDist(m_xt, cov_xt),
										ztDist(m_zt, cov_zt);

	for (int i = 0; i < nP; ++i) {
		points.push_back(Vector2d(xpDist(gen),zpDist(gen)));
	}

	for (int j = 0; j < mT; ++j) {

		Vector2d pose = Vector2d(xtDist(gen),ztDist(gen));
		std::vector<double> meas;

		for (int i = 0; i < nP; ++i)
			meas.push_back(project(points[i],pose));

		frames.push_back(frame{pose,meas});
	}

	std::vector<Vector2d> pGT = points;
	std::vector<frame> fGT = frames;

	gnOpt(points,frames,nIts);

	std::cout << "\n Poses: \n";
	for (int i = 0; i < mT; ++i){
		std::cout << fGT[i]._pose.transpose() << '\n';
		std::cout << frames[i]._pose.transpose() << "\n\n";
	}
	std::cout << "\n Points: \n";
	for (int i = 0; i < nP; ++i){
		std::cout << pGT[i].transpose() << '\n';
		std::cout << points[i].transpose() << "\n\n";
	}


}
