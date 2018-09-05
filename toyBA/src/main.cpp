/*
 * main.cpp
 *
 *  Created on: Sep 5, 2018
 *      Author: l_vis
 */

#include <vector>
#include <iostream>
#include <random>

#include <chrono>

#include "Eigen/Dense"

#include "ba.h"
#include "baWP.h"

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
					cov_zt = 0.1,
					cov_xtm = 0.1,
					cov_ztm = 0.1;
    std::default_random_engine gen;
    std::normal_distribution<double> 	xpDist(m_xp, cov_xp),
										zpDist(m_zp, cov_zp),
										xtDist(m_xt, cov_xt),
										ztDist(m_zt, cov_zt),
										xtmDist(0, cov_xtm),
										ztmDist(0, cov_ztm);

	for (int i = 0; i < nP; ++i) {
		points.push_back(Vector2d(xpDist(gen),zpDist(gen)));
	}


	for (int j = 0; j < mT; ++j) {

		Vector2d pose = Vector2d(xtDist(gen),ztDist(gen));
		std::vector<double> qMeas;
		std::vector<Vector2d> tMeas;

		for (int i = 0; i < nP; ++i)
			qMeas.push_back(project(points[i],pose));

		if(j == 0)
			tMeas.push_back(pose);// + Vector2d(xtmDist(gen),xtmDist(gen)));

		frames.push_back(frame{pose,tMeas,qMeas});
	}

	std::vector<Vector2d> pGT = points;
	std::vector<frame> fGT = frames;

	auto optPars = baWP::gnOpt(points,frames,nIts);

	std::cout << "\n Poses: \n";
	for (int i = 0; i < mT; ++i){
		std::cout << fGT[i]._t.transpose() << '\n';
		std::cout << frames[i]._t.transpose() << "\n\n";
	}
	std::cout << "\n Points: \n";
	for (int i = 0; i < nP; ++i){
		std::cout << pGT[i].transpose() << '\n';
		std::cout << points[i].transpose() << "\n\n";
	}

	std::cout << "nIts Error: " << optPars.first << '\n';
	std::cout << "Opt Error: " << optPars.second << '\n';

}
