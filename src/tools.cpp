#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    VectorXd rmse = VectorXd::Zero(4);
    if(estimations.size() <=0 or estimations.size() != ground_truth.size()){
        cout << "invalid estimation size"<<endl;
        return rmse;}

    else{
        
        for(int i=0; i < estimations.size(); ++i){
            VectorXd a =  estimations[i] - ground_truth[i];
            VectorXd b = a.array() * a.array();
            rmse = rmse + b;
            
            
        }
        rmse = rmse / estimations.size();
        rmse = rmse.array().sqrt();
        return rmse;
        
    }
}
