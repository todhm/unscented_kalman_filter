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
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;


  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ =  M_PI * 0.22;

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
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
  
  Nis_vec = VectorXd::Zero(2);
 
}

UKF::~UKF() {}


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
  MatrixXd A = P_aug.llt().matrixL();
  A = sqrt(lambda_ + n_aug_) * A;
  Xsig_aug.col(0) = x_aug;
  for(int i = 0 ; i < n_aug_; i ++){
    Xsig_aug.col(i + 1) = x_aug + A.col(i);
    Xsig_aug.col(i + n_aug_ + 1) = x_aug - A.col(i);

  }
  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, double delta_t){
  MatrixXd XsigPred = MatrixXd::Zero(n_x_,n_aug_*2 + 1);
  MatrixXd Xsig_aug = *Xsig_out;
  for(int i = 0; i< n_aug_ *2 +1 ; i ++){
    VectorXd predicted_state = Xsig_aug.col(i).head(n_x_);
    double division = Xsig_aug.col(i)(2) / Xsig_aug.col(i)(4) ;
    double sin_yawadd = sin(Xsig_aug.col(i)(3) + Xsig_aug.col(i)(4) * delta_t);
    double cos_yawadd = cos(Xsig_aug.col(i)(3) + Xsig_aug.col(i)(4) * delta_t);
    double sin_yaw = sin(Xsig_aug.col(i)(3));
    double cos_yaw = cos(Xsig_aug.col(i)(3));
    double nu_a = Xsig_aug.col(i)(5);
    double nu_yawadd = Xsig_aug.col(i)(6);
    
    predicted_state(0) += 0.5 * pow(delta_t, 2.0) * cos_yaw * nu_a;
    predicted_state(1) += 0.5 * pow(delta_t, 2.0) * sin_yaw * nu_a;
    predicted_state(2) += delta_t * nu_a;
    predicted_state(3) += 0.5 * pow(delta_t, 2.0)* nu_yawadd;
    predicted_state(4) += delta_t * nu_yawadd;
    
    if (fabs(Xsig_aug.col(i)(4)) <0.001 ){
      predicted_state(0) += Xsig_aug.col(i)(2) * cos_yaw * delta_t  ;
      predicted_state(1) += Xsig_aug.col(i)(2) * sin_yaw * delta_t;
    }
    else{
      predicted_state(0) += division * (sin_yawadd - sin_yaw);
      predicted_state(1) += division * ( -cos_yawadd + cos_yaw);
      predicted_state(3) += Xsig_aug.col(i)(4)* delta_t;
    }
    XsigPred.col(i) = predicted_state;
  }
  Xsig_pred_ = XsigPred;

};

void UKF::PredictMeanAndCovariance(){
  
  MatrixXd Ppred = MatrixXd::Zero(n_x_, n_x_);
  
  //weights_.fill(0.5/(n_aug_+lambda_));
  //weights_(0) = lambda_ / (lambda_ + n_aug_);
  
  x_ = Xsig_pred_* weights_;
  
  for(int i =0; i < 2 * n_aug_ +1 ; i ++){
    VectorXd xdiff =Xsig_pred_.col(i) - x_;
    xdiff[3] -=  (2* M_PI) * floor((xdiff[3] + M_PI) / (2*M_PI));
    Ppred += weights_(i) * xdiff * xdiff.transpose();
  }
  P_ = Ppred;
  
}
void UKF::Prediction(double delta_t){
  MatrixXd Xsig = MatrixXd::Zero(n_aug_, n_aug_*2 + 1);
  UKF::AugmentedSigmaPoints(&Xsig);
  UKF::SigmaPointPrediction(&Xsig ,delta_t);
  UKF::PredictMeanAndCovariance();
  
 }


void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* Zsig, MatrixXd* S_out){
  MatrixXd R = MatrixXd::Zero(3,3);
  VectorXd zk_1 = VectorXd(3) ;
  MatrixXd S_1 = MatrixXd::Zero(3,3);
  R(0,0) = std_radr_ * std_radr_;
  R(1,1) = std_radphi_ * std_radphi_;
  R(2,2) = std_radrd_ * std_radrd_;
  MatrixXd Zk_1 = MatrixXd::Zero(3,n_aug_*2 + 1 );
  for(int i = 0 ; i< n_aug_ * 2 + 1; i ++ ){
    float px  = Xsig_pred_.col(i)(0);
    float py  = Xsig_pred_.col(i)(1);
    float v = Xsig_pred_.col(i)(2);
    float yaw = Xsig_pred_.col(i)(3);
    if(fabs(px) < 0.001){
      px = 0.001;
    }
    if(fabs(py) < 0.001){
      py = 0.001;
    }
    
    float rho = sqrt(px * px + py * py);
    if(rho < 0.001){
      rho = 0.001;
    }

    float theta = atan2(py, px);

    float rho_dot = (px * cos(yaw) *v  + py * sin(yaw) * v) / rho ;
    
    Zk_1.col(i) << rho, theta, rho_dot;
  }
  zk_1 = Zk_1 * weights_;
  for(int i = 0 ; i < n_aug_ * 2 +1 ; i ++){
    VectorXd zdiff =Zk_1.col(i) - zk_1;
    zdiff[1] -=  (2* M_PI) *floor( (zdiff[1] + M_PI) /(2*M_PI));
   
    S_1 += weights_(i) * zdiff * zdiff.transpose() ;
  
  }
  S_1 += R;
 
  *z_out = zk_1;
  *S_out = S_1;
  *Zsig = Zk_1; 
  
}
void UKF::PredictLidarMeasurement(VectorXd* z_out,MatrixXd* Zsig ,MatrixXd* S_out){
  MatrixXd R = MatrixXd::Zero(2,2);
  VectorXd zk_1;
  MatrixXd S_1 = MatrixXd::Zero(2,2);
  R(0,0) = std_laspx_;
  R(1,1) = std_laspy_ ;
  MatrixXd Zk_1 = Xsig_pred_.topLeftCorner(2,n_aug_*2 + 1);
  zk_1 = Zk_1 * weights_;
  for(int i = 0 ; i < n_aug_ * 2 +1 ; i ++){
    VectorXd zdiff =Zk_1.col(i) - zk_1;
    S_1 += weights_(i) * zdiff * zdiff.transpose() ;
    
  }
  S_1 += R;

  *z_out = zk_1;
  *S_out = S_1;
  *Zsig = Zk_1;
  

}


void UKF::UpdateState(VectorXd* z_out,MatrixXd* Zsig, MatrixXd* S_out,MeasurementPackage meas_package){
  VectorXd measurement = meas_package.raw_measurements_;
  MatrixXd T;
  if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    T = MatrixXd::Zero(n_x_, n_x_ -3);
  }
  else{
    T = MatrixXd::Zero(n_x_,n_x_- 2);
  }
  MatrixXd zpred = *z_out;
  MatrixXd Spred = *S_out;
  MatrixXd Zsig_ = *Zsig;

  for(int i = 0; i< 2*n_aug_ + 1 ; i++){
    VectorXd zdiff = Zsig_.col(i) - zpred;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    zdiff[1] -=  (2* M_PI) *floor( (zdiff[1] + M_PI) /(2*M_PI));
    
    }
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    xdiff(3) -=  (2* M_PI) *floor( (xdiff(3) + M_PI) /(2*M_PI));
    //cout<<xdiff * zdiff.transpose()<<endl;
    T += weights_(i) * (xdiff * zdiff.transpose());

  }
  MatrixXd K ;
  K = T * Spred.inverse();
  VectorXd zdiff_direct =measurement - zpred;
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
  zdiff_direct(1) -= 2 *M_PI * floor((zdiff_direct(1) + M_PI) / (2 * M_PI));
  }
  x_ = x_ + K*(zdiff_direct);
  P_ = P_ - K*Spred*K.transpose();
  double NIS;
  NIS = zdiff_direct.transpose() * Spred.inverse() * zdiff_direct;
  //add total number of NIS HAVE been calculated so far.
  Nis_vec(0) += 1;
  //find a case where NIS is beyond 95% distribution range.
  if(NIS > 7.8){
    Nis_vec(1) += 1;
  }
  cout<<"NIS: "<<NIS<<endl;
  cout<<"NIS Percent: " <<Nis_vec(1) / Nis_vec(0)<<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  VectorXd z;
  MatrixXd S;
  MatrixXd Zsig;
  VectorXd measurement = meas_package.raw_measurements_;
  UKF::PredictRadarMeasurement(&z,&Zsig, &S);
  UKF::UpdateState(&z, &Zsig, &S, meas_package);
  
  
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
  MatrixXd Zsig;
  VectorXd measurement = meas_package.raw_measurements_;
  UKF::PredictLidarMeasurement(&z,&Zsig, &S);
  UKF::UpdateState(&z,&Zsig,&S, meas_package );
  
}





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
      
      float rho = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      float rhodot = meas_package.raw_measurements_(2);
      

      float px = rho * cos(theta);
      float py = rho * sin(theta);
      float vx = rhodot * cos(theta);
      float vy = rhodot * sin(theta);

      x_ <<  px,py,sqrt(vx*vx +vy*vy),0,0;
    }
    is_initialized_ = true;
    return ;
    
  }
  
  else{
    double delta_t = (meas_package.timestamp_- time_us_ )/1000000.0;;
    time_us_ = meas_package.timestamp_;
    while (delta_t > 0.1)
    {
      const double dt = 0.05;
      Prediction(dt);
      delta_t -= dt;
    }
    UKF::Prediction(delta_t);
    

    
    if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
      UKF::UpdateLidar(meas_package);

    }
    else if ( meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ ) {
      if(use_radar_ == true){
      UKF::UpdateRadar(meas_package);
      }

    }
    time_us_ = meas_package.timestamp_;
    return;
    
  }
  
}

