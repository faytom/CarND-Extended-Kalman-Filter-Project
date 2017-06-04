#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd RMSE(4);

	VectorXd estLast(4);
	estLast = estimations.back();
	VectorXd gtLast(4);
	gtLast = ground_truth.back();

	RMSE(0) = sqrt((estLast(0) - gtLast(0)) *  (estLast(0) - gtLast(0)));
	RMSE(1) = sqrt((estLast(1) - gtLast(1)) *  (estLast(1) - gtLast(1)));
	RMSE(2) = sqrt((estLast(2) - gtLast(2)) *  (estLast(2) - gtLast(2)));
	RMSE(3) = sqrt((estLast(3) - gtLast(3)) *  (estLast(3) - gtLast(3)));
	return RMSE;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
