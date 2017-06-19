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
    //initialization
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;
    
    //check if the size of data sets matches
    if ((estimations.size() != ground_truth.size()) || estimations.size() == 0) {
        cout << "Error in estimations or ground truth data" << endl;
        return rmse;
    }
    
    //accumulate the squared residuals
    for (int i=0; i < estimations.size(); ++i) {
        VectorXd residual = estimations[i] - ground_truth[i];
        
        //element wise multiplication
        residual = residual.array() * residual.array();
        rmse += residual;
    }
    
    //calculate the mean
    rmse = rmse / estimations.size();
    
    //calculate the square root
    rmse = sqrt(rmse.array());
    
    //return the result
    return rmse;
}
