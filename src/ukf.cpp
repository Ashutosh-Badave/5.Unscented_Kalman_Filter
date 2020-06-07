#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include<fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 4;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 8;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;

  Xsig_pred_ = MatrixXd(5,15);

  n_x_ = 5;

  n_aug_ = 7;
  
  lambda_ = 3 - n_aug_;
  
  weights_ = VectorXd(15);
       
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_)
  {
    P_ = MatrixXd::Identity(5,5);
    if(meas_package.sensor_type_ == MeasurementPackage::LASER) 
    {
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0,
            0,
            0;

      // Init Process covariance matrix 
      P_(0,0) = std_laspx_*std_laspx_;
      P_(1,1) = std_laspy_*std_laspy_;
      
    }
    else 
    {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rhodot = meas_package.raw_measurements_[2];

      double px = rho * cos(phi);
      double py = rho * sin(phi);
      x_ << px,
            py,
            rhodot,
            phi,
            0;
      // Init Process covariance matrix 
      P_(0,0) = std_radr_*std_radr_;
      P_(1,1) = std_radr_*std_radr_;
      P_(2,2) = std_radr_*std_radr_;
    
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0L;
  time_us_ = meas_package.timestamp_;

  //Prediction of next position
  Prediction(dt);
  
  //Update depending on sensor type
  if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }
  
}

void UKF::Prediction(double delta_t) {


//--------------------1. Generate Augmented sigma points--------------------------
 
  MatrixXd Xsig_Aug = MatrixXd(n_aug_,2*n_aug_+1);
  
  // Augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // augmented state variable
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);

  // Augmented mean state
  x_aug.head(5) = x_;
  x_aug[5] = 0; 
  x_aug[6] = 0;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = pow(std_a_,2);
  P_aug(6,6) = pow(std_yawdd_,2);
  
  // create augmented sigma points
  MatrixXd A = P_aug.llt().matrixL();

  Xsig_Aug.fill(0.0);
  // first coloum of sigma points
  Xsig_Aug.col(0) = x_aug;

  for (int i = 0; i<n_aug_;i++)
  {
    Xsig_Aug.col(i+1) = x_aug + sqrt(n_aug_ + lambda_) * A.col(i);
    Xsig_Aug.col(i+1+n_aug_) = x_aug - sqrt(n_aug_ + lambda_) * A.col(i);
  }
  
//--------------------2. Predict Sigma Point------------------------------
  // predict sigma points
  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_+1);
  Xsig_pred_.fill(0.0);
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    // extract values for better readability
    double p_x = Xsig_Aug(0,i);
    double p_y = Xsig_Aug(1,i);
    double v = Xsig_Aug(2,i);
    double yaw = Xsig_Aug(3,i);
    double yawd = Xsig_Aug(4,i);
    double nu_a = Xsig_Aug(5,i);
    double nu_yawdd = Xsig_Aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
//--------------------3. Predict the mean and covariance----------------------------------
 // vector for weights
  weights_ = VectorXd::Zero(2*n_aug_+1);
  
 // set weights
  for(int i=0;i<2*n_aug_+1;i++)
  {
      if(i==0)
      {
          weights_(i)=lambda_/(lambda_+n_aug_);
          
      }
      else
      {
          weights_(i)=(0.5/(lambda_+n_aug_));
      }
  }
  
 // predicted state mean
 x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    x_ +=  weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  // print result
 // std::cout << "Predicted state" << std::endl;
  //std::cout << x_ << std::endl;
  //std::cout << "Predicted covariance matrix" << std::endl;
 // std::cout << P_ << std::endl;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {

//----------------Convert the prediction to measurement space--------------------
  int n_z = 2;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  

  VectorXd z_pred = VectorXd(n_z);
  
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  MatrixXd R = MatrixXd(n_z,n_z);
  R(0,0) = pow(std_laspx_,2);
  R(1,1) = pow(std_laspy_,2);

  Zsig.fill(0.0);
  for (int i=0;i<2 * n_aug_+1;++i)
  {
      Zsig(0,i) = Xsig_pred_(0,i);
      Zsig(1,i) = Xsig_pred_(1,i);
  }
  
//----------------1. Predict Measurement-----------------------------
  z_pred.fill(0.0);
// calculate mean predicted measurement
  for (int i=0;i<2 * n_aug_+1;++i)
  {
      z_pred += weights_(i)*Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  S.fill(0.0);
  Tc.fill(0.0);
  for (int i=0;i<2 * n_aug_+1;i++)
  {
      VectorXd z_diff = Zsig.col(i) - z_pred;

      S += (weights_(i)* z_diff * z_diff.transpose());
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  S += R;
  // print result
 // std::cout << "z_pred: " << std::endl << z_pred << std::endl;
 // std::cout << "S: " << std::endl << S << std::endl;

//----------------2. Update -------------------------------------
  // read the measurement data
  VectorXd z = meas_package.raw_measurements_;
  
  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  x_ += K*(z-z_pred);
  P_ -= K*S*K.transpose();

  // print result
  //std::cout << "NIS value: " << std::endl << (z-z_pred).transpose() * S.inverse() * (z-z_pred) << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  std::ofstream log("../NIS_Laser.txt",std::ios_base::app | std::ios_base::out);
  log<< (z-z_pred).transpose() * S.inverse() * (z-z_pred);
  log << "\n";

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {

//------------------Convert to measurement space---------------------
  int n_z = 3;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0.0);
  R(0,0) = pow(std_radr_,2);
  R(1,1) = pow(std_radphi_,2);
  R(2,2) = pow(std_radrd_,2);
  // transform sigma points into measurement space
  Zsig.fill(0.0);
  for (int i=0;i<2 * n_aug_+1;i++)
  {
      double pos_x = Xsig_pred_(0,i);
      double pos_y = Xsig_pred_(1,i);
      double vel = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);
      
      double rad = sqrt(pow(pos_x,2)+pow(pos_y,2));
      double ang = atan2(pos_y,pos_x);
      double radd = (pos_x*cos(yaw)*vel + pos_y*sin(yaw)*vel) / rad;
      
      Zsig(0,i) = rad;
      Zsig(1,i) = ang;
      Zsig(2,i) = radd;
  }

//----------------1. Predict measurement ---------------------------  
// calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0;i<2 * n_aug_+1;i++)
  {
      z_pred += weights_(i)*Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  S.fill(0.0);
  Tc.fill(0.0);
  for (int i=0;i<2 * n_aug_+1;i++)
  {
      //measurement difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    //state difference 
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    S += (weights_(i)* z_diff * z_diff.transpose());

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  S += R;

//----------------2. Update ------------------------------
  // read the measurement data
  VectorXd z = meas_package.raw_measurements_;
  /*z << meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1),
       meas_package.raw_measurements_(2); 
*/
  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z-z_pred;

  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ += K*z_diff;
  P_ -= K*S*K.transpose();

  // print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
 // std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
 //std::cout << "NIS value Radar: " << std::endl << (z-z_pred).transpose() * S.inverse() * (z-z_pred) << std::endl;
  std::ofstream log_R("../NIS_Radar.txt",std::ios_base::app | std::ios_base::out);
  log_R<< (z-z_pred).transpose() * S.inverse() * (z-z_pred);
  log_R << "\n";
}