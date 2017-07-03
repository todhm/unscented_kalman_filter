#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;




/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF(){
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;


  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ =  M_PI / 8;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  is_initialized_ = false;
  
    
  //Set time to 0
  time_us_ = 0 ;
    
   ///* State dimension
  n_x_ = 5;
  
  // initial state vector
  x_ = VectorXd::Zero(n_x_);
  
  // initial covariance matrix
  P_ = MatrixXd::Identity(5,5);

    
  //Augmented state dimension
  n_aug_  = 7;
  lambda_ = 3 - n_x_ ;
  
  Xsig_pred_ =  MatrixXd::Zero(n_x_,2*n_aug_ + 1);
  
  
  weights_ = VectorXd::Zero(2*n_aug_ + 1);
  weights_.fill(1/(2 *(lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  
}
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */


UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    if( is_initialized_ != true){
        time_us_ = meas_package.timestamp_;

      
     
        cout<<"Initialized"<<endl;
      
        if (meas_package.sensor_type_ == MeasurementPackage::LASER){
            float px =  meas_package.raw_measurements_(0);
            float py = meas_package.raw_measurements_(1);
            x_(0) = px;
            x_(1) = py;
    
          
    
        }
        else{

            float ro = meas_package.raw_measurements_(0);
            float theta = meas_package.raw_measurements_(1);
            float px = ro * cos(theta);
            float py = ro * sin(theta);
            x_ << px,py,0,0,0;
          
        }

    
      return ;
    }
    else{
      double delta_t =meas_package.timestamp_ - time_us_;
      UKF::Prediction(delta_t);
      time_us_ = meas_package.timestamp_;
      
    
    
      if(meas_package.sensor_type_ == MeasurementPackage::LASER){
        UKF::UpdateLidar(meas_package);
        
      }
      else{
        UKF::UpdateRadar(meas_package);
      }

    }
    /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */




void UKF::AugmentedSigmaPoints(MatrixXd *Xsig_out){
  
  
  MatrixXd Xsig_aug =MatrixXd::Zero(n_aug_, 2 *n_aug_ + 1);
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd Q_ = MatrixXd(2,2);
  
  Q_.fill(0);
  Q_(0,0) = std_a_ * std_a_;
  Q_(1,1) = std_yawdd_ * std_yawdd_;

  x_aug.head(n_x_) = x_;
  x_aug.tail(n_aug_ - n_x_) << 0,0;
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = Q_;
  MatrixXd A = P_.llt().matrixL();
  A = sqrt(lambda_ + n_aug_) * A;
  Xsig_aug.col(0) = x_aug;
  for(int i = 0 ; i < n_aug_; i ++){
    Xsig_aug.col(i + 1) = x_ + A.col(i);
    Xsig_aug.col(i +1)(3) -= 2*M_PI* floor((Xsig_aug.col(i +1)(3) + M_PI)/ (2*M_PI));
    Xsig_aug.col(i + n_aug_ + 1) = x_aug - A.col(i);
    Xsig_aug.col(i + n_aug_ + 1)(3) -= 2*M_PI* floor((Xsig_aug.col(i + n_aug_ +1)(3) + M_PI)/ (2*M_PI));

  }
  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, double delta_t){
  MatrixXd XsigPred = MatrixXd::Zero(n_x_,n_aug_*2 + 1);
  MatrixXd Xsig_aug = *Xsig_out;
  for(int i = 0; i< n_aug_ ; i ++){
    VectorXd predicted_state = Xsig_aug.col(i).head(n_x_);
    double division = predicted_state(2) / predicted_state(4) ;
    double sin_yawadd = sin(predicted_state(3) + predicted_state(4) * delta_t);
    double cos_yawadd = cos(predicted_state(3) + predicted_state(4) * delta_t);
    double sin_yaw = sin(predicted_state(3));
    double cos_yaw = cos(predicted_state(3));
    double noise_acc = Xsig_aug.col(i)(5);
    double noise_yawdot = Xsig_aug.col(i)(6);
    
    predicted_state(0) += 0.5 * pow(delta_t, 2.0) * cos_yaw * noise_acc;
    predicted_state(1) += 0.5 * pow(delta_t, 2.0) * sin_yaw * noise_acc;
    predicted_state(2) += delta_t * noise_acc;
    predicted_state(3) += 0.5 * pow(delta_t, 2.0)* noise_yawdot;
    predicted_state(4) += delta_t * noise_yawdot;
    
    if (predicted_state(4) ==0 ){
      predicted_state(0) += predicted_state(2) * cos_yaw * delta_t ;
      predicted_state(1) += predicted_state(2) * sin_yaw * delta_t;
    }
    else{
      predicted_state(0) += division * (sin_yawadd - sin_yaw);
      predicted_state(1) += division * ( -cos_yawadd + cos_yaw);
      predicted_state(3) += predicted_state(4) * delta_t;
    }
    XsigPred.col(i) = predicted_state;
  
  }
  Xsig_pred_ = XsigPred;
};

void UKF::PredictMeanAndCovariance(){
  
  MatrixXd Ppred = MatrixXd::Zero(n_x_, n_x_);

  
  x_ = Xsig_pred_* weights_;
  for(int i =0; i < 2 * n_aug_ +1 ; i ++){
    Ppred += weights_(i) * (Xsig_pred_.col(i) - x_) * (Xsig_pred_.col(i) - x_).transpose();
  }
  P_ = Ppred;
  
}
void UKF::Prediction(double delta_t){
  MatrixXd Xsig = MatrixXd::Zero(n_aug_, n_aug_*2 + 1);
  cout<<"predict"<<endl;
  UKF::AugmentedSigmaPoints(&Xsig);
  UKF::SigmaPointPrediction(&Xsig ,delta_t);
  UKF::PredictMeanAndCovariance();
  
 }


void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out){
  MatrixXd R = MatrixXd::Zero(3,3);
  VectorXd zk_1 ;
  MatrixXd S_1 = MatrixXd::Zero(3,3);
  R(0,0) = std_radr_ * std_radr_;
  R(1,1) = std_radphi_ * std_radphi_;
  R(2,2) = std_radrd_ * std_radrd_;
  MatrixXd Zk_1 = MatrixXd::Zero(3,n_aug_*2 + 1 );
  for(int i = 0 ; i< n_aug_ * 2 + 1; i ++ ){
    float px  = Xsig_pred_(0);
    float py  = Xsig_pred_(1);
    float v = Xsig_pred_(2);
    float yaw = Xsig_pred_(3);
    float rho = sqrt(px * px + py * py);
    float theta = atan2(py, px);
    float rho_dot = (px * cos(yaw) *v  + py * sin(yaw) * v) / rho ;
    Zk_1.col(i) << rho, theta, rho_dot;
  }
   zk_1 = Zk_1 * weights_;
  for(int i = 0 ; i < n_aug_ * 2 +1 ; i ++){
    VectorXd zdiff =Zk_1.col(i) - zk_1;
    zdiff(1) -=  2* M_PI *floor( (zdiff(1) + M_PI) /(2*M_PI));
    S_1 += weights_(i) * zdiff * zdiff.transpose() + R;
  
  }
  *z_out = zk_1;
  *S_out = S_1;
  
  
}
void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out){
  MatrixXd R = MatrixXd::Zero(2,2);
  VectorXd zk_1;
  MatrixXd S_1 = MatrixXd::Zero(2,2);
  R(0,0) = std_laspx_;
  R(1,1) = std_laspy_ ;
  MatrixXd Zk_1 = Xsig_pred_.topLeftCorner(2,n_aug_*2 + 1);
  zk_1 = Zk_1 * weights_;
  for(int i = 0 ; i < n_aug_ * 2 +1 ; i ++){
    VectorXd zdiff =Zk_1.col(i) - zk_1;
    S_1 += weights_(i) * zdiff * zdiff.transpose() + R;
    
  }
  *z_out = zk_1;
  *S_out = S_1;

}


void UKF::UpdateState(VectorXd* z_out, MatrixXd* S_out,MeasurementPackage meas_package){
  VectorXd measurement = meas_package.raw_measurements_;
  MatrixXd T = MatrixXd::Zero(n_x_,n_x_ -2);
  MatrixXd zpred = *z_out;
  MatrixXd Spred = *S_out;
  for(int i = 0; i< 2*n_aug_ + 1 ; i++){
    VectorXd zdiff = zpred.col(i)  - measurement;
    if (meas_package.sensor_type_ == MeasurementPackage::LASER){
    zdiff(1) -=  2* M_PI *floor( (zdiff(1) + M_PI) /(2*M_PI));
    }
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    T += weights_(i) * xdiff * zdiff.transpose();

  }
  MatrixXd K ;
  K = T*Spred.inverse();
  x_ = x_ + K*(measurement - zpred);
  P_ = P_ - K*Spred*K.transpose();
  
  
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  VectorXd z;
  MatrixXd S;
  VectorXd measurement = meas_package.raw_measurements_;
  UKF::PredictRadarMeasurement(&z, &S);
  UKF::UpdateState(&z,&S, meas_package);
  
  
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  VectorXd z;
  MatrixXd S;
  VectorXd measurement = meas_package.raw_measurements_;
  UKF::PredictRadarMeasurement(&z, &S);
  UKF::UpdateState(&z,&S, meas_package);
  
  /**
   TODO:
   
   Complete this function! Use lidar data to update the belief about the object's
   position. Modify the state vector, x_, and covariance, P_.
   
   You'll also need to calculate the lidar NIS.
   */
  
}
