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
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
    // State dimension
    n_x_ = 5;
    
    // Augmented state dimension
    n_aug_ = 7;
    
    // Sigma point spreading parameter
    lambda_ = 3 - n_aug_;
    
    // is initialized or not?
    is_initialized_ = false;
    
    //time
    previous_timestamp_ = 0.0;
    
    // predicted augmented matrix
    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
    
    // weights of sigma points
    weights_ = VectorXd(2*n_aug_ + 1);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i=1; i < 2*n_aug_ + 1; i++) {
        weights_(i) = 0.5 / (lambda_ + n_aug_);
    }
    
    // initialize NIS for Lidar
    NIS_laser_ = 0.0;
    
    // initialize NIS for radar
    NIS_radar_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    if (!is_initialized_) {
        
        // covariance
        P_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1;
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            // initialize state vector with first measurement data from radar
            float px = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
            float py = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
            x_ << px, py, 0, 0, 0;
            previous_timestamp_ = meas_package.timestamp_;
        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            // initialize state vector with first measurement data from lidar
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
            previous_timestamp_ = meas_package.timestamp_;
    }
    
    // Done initialization
    is_initialized_ = true;
    return;
}


    float delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = meas_package.timestamp_;
    
    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        Prediction(delta_t);
        UpdateRadar(meas_package);
    } else if(use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
        Prediction(delta_t);
        UpdateLidar(meas_package);
    }
}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    //augment 2 acceleration with state vector
    VectorXd x_aug_ = VectorXd(n_aug_);
    x_aug_.head(n_x_) = x_;
    x_aug_(5) = 0.0;
    x_aug_(6) = 0.0;
    
    //augment 2 noise covariance with the state covariance matrix
    MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(n_x_, n_x_) = P_;
    P_aug_(5, 5) = std_a_ * std_a_;
    P_aug_(6, 6) = std_yawdd_ * std_yawdd_;
    
    //square root of P_aug_
    MatrixXd P_aug_sqrt_ = P_aug_.llt().matrixL();
    
    //generate sigma points around current state
    MatrixXd Xsig_ = MatrixXd(n_aug_, 2*n_aug_ + 1);
    Xsig_.col(0) = x_aug_;
    for (int i=0; i < n_aug_; i++) {
        Xsig_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_) * P_aug_sqrt_.col(i);
        Xsig_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * P_aug_sqrt_.col(i);
    }
    
    //sigma points prediction
    for (int i=0; i < 2*n_aug_ + 1; i++) {
        
        //unpack the current augmented state vector
        double p_x_ = Xsig_(0, i);
        double p_y_ = Xsig_(1, i);
        double v_ = Xsig_(2, i);
        double yaw_ = Xsig_(3, i);
        double yaw_rate_ = Xsig_(4, i);
        double acce_noise = Xsig_(5, i);
        double yaw_noise = Xsig_(6, i);
        
        //calculate noise vector
        VectorXd noise_ = VectorXd(5);
        double delta_t_sq_ = delta_t * delta_t;
        
        noise_(0) = delta_t_sq_*cos(yaw_)*acce_noise/2;
        noise_(1) = delta_t_sq_*sin(yaw_)*acce_noise/2;
        noise_(2) = delta_t*acce_noise;
        noise_(3) = delta_t_sq_*yaw_noise/2;
        noise_(4) = delta_t*yaw_noise;
        
        if (fabs(yaw_rate_) > 0.001) {
            Xsig_pred_(0, i) = p_x_ + (sin(yaw_+yaw_rate_*delta_t)-sin(yaw_))*v_/yaw_rate_ + noise_(0);
            Xsig_pred_(1, i) = p_y_ + (-cos(yaw_+yaw_rate_*delta_t) + cos(yaw_))*v_/yaw_rate_ + noise_(1);
        } else {
            Xsig_pred_(0, i) = p_x_ + v_*cos(yaw_)*delta_t + noise_(0);
            Xsig_pred_(1, i) = p_y_ + v_*sin(yaw_)*delta_t + noise_(1);
        }
        
        Xsig_pred_(2, i) = v_ + noise_(2);
        Xsig_pred_(3, i) = yaw_ + yaw_rate_*delta_t + noise_(3);
        Xsig_pred_(4, i) = yaw_rate_ + noise_(4);
    }
    
    // Predict mean and covariance
    x_.fill(0.0);//***************************************************
    for (int i=0; i < 2*n_aug_ + 1; i++) {
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }
    
    VectorXd diff_to_mean = VectorXd(n_x_);
    
    P_.fill(0.0);
    for (int i=0; i < 2*n_aug_ + 1; i++) {
        diff_to_mean = Xsig_pred_.col(i) - x_;
        
        // check if yaw is between -pi and pi
        while (diff_to_mean(3) > M_PI) {
            diff_to_mean(3) -= 2*M_PI;
        }
        while (diff_to_mean(3) < -M_PI) {
            diff_to_mean(3) += 2*M_PI;
        }
        
        P_ = P_ + weights_(i) * diff_to_mean * diff_to_mean.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    //Take in the real measurement
    int n_w_ = 2;
    VectorXd w_ = VectorXd(n_w_);
    w_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
    
    //measurement covariance matrix
    MatrixXd R_ = MatrixXd(n_w_, n_w_);
    R_ << std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
    
    //H
    MatrixXd H_ = MatrixXd(n_w_, n_x_);
    H_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0;
    
    //measurement update
    VectorXd Wpred_ = H_ * x_;
    VectorXd y = w_ - Wpred_;
    MatrixXd Ht = H_.transpose();
    MatrixXd PHt = P_ * Ht;
    MatrixXd S_ = H_* PHt + R_;
    MatrixXd Si = S_.inverse();
    MatrixXd K_ = PHt * Si;
    
    x_ = x_ + (K_ * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K_*H_) * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    int n_z_ = 3;
    VectorXd z_ = VectorXd(n_z_);
    z_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
    
    // transfer all predicted sigma point to polar coordinates
    MatrixXd Zsig_ = MatrixXd(n_z_, 2*n_aug_ + 1);
    
    for (int i=0; i < 2*n_aug_ + 1; i++) {
        
        double p_x_ = Xsig_pred_(0, i);
        double p_y_ = Xsig_pred_(1, i);
        double v_ = Xsig_pred_(2, i);
        double yaw_ = Xsig_pred_(3, i);
        double yaw_rate_ = Xsig_pred_(4, i);
        
        if (fabs(p_x_) < 0.001) {
            p_x_ = 0.001;
        }
        
        double rho = sqrt(p_x_*p_x_ + p_y_*p_y_);
        if (fabs(rho) < 0.0001) {
            rho = 0.0001;
        }
        
        Zsig_(0, i) = rho;
        Zsig_(1, i) = atan2(p_y_, p_x_);
        Zsig_(2, i) = (p_x_*cos(yaw_)+ p_y_*sin(yaw_))*v_ / rho;
    }
    
    //Calculate Predicted Measurement Mean
    VectorXd Zpred_ = VectorXd(n_z_);
    for (int i=0; i < 2*n_aug_ + 1; i++) {
        Zpred_ = weights_(i) * Zsig_.col(i);
    }
    
    //Measurement covariance
    MatrixXd R_ = MatrixXd(n_z_, n_z_);
    R_ << std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0, std_radrd_*std_radrd_;
    
    //Calculate Predicted Measurement Covariance
    MatrixXd S_ = MatrixXd(n_z_, n_z_);
    
    S_.fill(0.0);
    VectorXd diff_to_mean = VectorXd(n_z_);
    
    for (int i=0; i < 2*n_aug_ + 1; i++) {
        diff_to_mean = Zsig_.col(i) - Zpred_;
        
        // check if phi is between -pi and pi
        while (diff_to_mean(1) > M_PI) {
            diff_to_mean(1) -= 2*M_PI;
        }
        while (diff_to_mean(1) < -M_PI) {
            diff_to_mean(1) += 2*M_PI;
        }
        
        S_ = S_ + weights_(i) * diff_to_mean * diff_to_mean.transpose();
        
    }
    
    //add measurement covariance matrix
    S_ = S_ + R_;
    
    //cross-correlation matrix
    VectorXd x_diff_to_mean = VectorXd(n_x_);
    MatrixXd T_ = MatrixXd(n_x_, n_z_);
    T_.fill(0.0);
    
    for (int i=0; i < 2*n_aug_ + 1; i++) {
        x_diff_to_mean = Xsig_pred_.col(i) - x_;
        
        // check if yaw is between -pi and pi
        while (x_diff_to_mean(3) > M_PI) {
            x_diff_to_mean(3) -= 2*M_PI;
        }
        while (x_diff_to_mean(3) < -M_PI) {
            x_diff_to_mean(3) += 2*M_PI;
        }
        
        T_ = T_ + weights_(i) * x_diff_to_mean * diff_to_mean.transpose();
    }
    
    //calculate Kalman gain
    MatrixXd K_ = T_ * S_.inverse();
    
    //calculate NIS value
    NIS_radar_ = (z_ - Zpred_).transpose() * S_.inverse() * (z_ - Zpred_);
    
    //update state and covariance matrix
    x_ = x_ + K_*(z_ - Zpred_);
    P_ = P_ - K_*S_*K_.transpose();
}
