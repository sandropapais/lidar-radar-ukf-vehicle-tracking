#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5,5);

  // time when the state is true, in us
  time_us_ = 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
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
  
  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i = 1; i < 2*n_aug_+1; ++i) {
    weights_(i) = 0.5/(n_aug_+lambda_);
  }

  // Predicted sigma point matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2*n_aug_+1);

  // Debug flag
  debug_flg = false;

  // Print results flag
  print_flg = true;
  if(print_flg){
		FILE *fpt;
		fpt = fopen("nis_out.csv", "w");
	  fprintf(fpt, "sen_ID_lid_rad, time_s, nis\n");
		fclose(fpt);
 		fpt = fopen("mes_out.csv", "w");
	  fprintf(fpt, "sen_ID_lid_rad, time_s, x, y, v\n");
		fclose(fpt);
	}

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  
  if(debug_flg) std::cout << ">>>>> START OF MEASURMENT PROCESSING <<<<<" << std::endl;
  if(debug_flg) std::cout << "timestamp (s) = " << meas_package.timestamp_/1000000.0 << std::endl;  
  // Initialize state with first measurement
  if (is_initialized_ == 0)
  {
    if(debug_flg) std::cout << "First measurement" << std::endl;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      if(debug_flg) std::cout << "Lidar measurement: " << meas_package.raw_measurements_ << std::endl;
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0, std_laspy_*std_laspy_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
      if(print_flg){
 	  	  FILE *fpt;
 		    fpt = fopen("mes_out.csv", "a");
	      fprintf(fpt, "%d, %f, %f, %f, %f\n", 0, meas_package.timestamp_/1e6, meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0.0);
    		fclose(fpt);
      }
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      if(debug_flg) std::cout << "Radar measurement: " << meas_package.raw_measurements_ << std::endl;
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_d = meas_package.raw_measurements_(2);
      x_ << rho*cos(phi), rho*sin(phi), rho_d, 0, 0;
      P_ << std_radr_*std_radr_, 0, 0, 0, 0,
            0, std_radr_*std_radr_, 0, 0, 0,
            0, 0, std_radrd_*std_radrd_, 0, 0,
            0, 0, 0, std_radphi_*std_radphi_, 0,
            0, 0, 0, 0, 1;
      if(print_flg){
 	  	  FILE *fpt;
 		    fpt = fopen("mes_out.csv", "a");
	      fprintf(fpt, "%d, %f, %f, %f, %f\n", 1, meas_package.timestamp_/1e6, rho*cos(phi), rho*sin(phi), rho_d);
    		fclose(fpt);
      }

    }
    
    if(debug_flg) std::cout << "State vector initialized:" << x_ << std::endl;

    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }
  
  // Predict step
  double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);

  // Lidar correct step
  if (meas_package.sensor_type_ == MeasurementPackage::LASER &&
      use_laser_)
  {
    if(print_flg){
 	    FILE *fpt;
 		  fpt = fopen("mes_out.csv", "a");
	    fprintf(fpt, "%d, %f, %f, %f, %f\n", 0, meas_package.timestamp_/1e6, meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0.0);
    	fclose(fpt);
    }
    if(debug_flg) std::cout << "Lidar measurement: " << meas_package.raw_measurements_ << std::endl;
    UpdateLidar(meas_package);
    if(debug_flg) std::cout << "State vector corrected:" << x_ << std::endl;
  }
  // Radar correct step
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR &&
      use_radar_)
  {
    if(print_flg){
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_d = meas_package.raw_measurements_(2);
  	  FILE *fpt;
	    fpt = fopen("mes_out.csv", "a");
      fprintf(fpt, "%d, %f, %f, %f, %f\n", 1, meas_package.timestamp_/1e6, rho*cos(phi), rho*sin(phi), rho_d);
  		fclose(fpt);
    }
    if(debug_flg) std::cout << "Radar measurement: " << meas_package.raw_measurements_ << std::endl;
    UpdateRadar(meas_package);
    if(debug_flg) std::cout << "State vector corrected:" << x_ << std::endl;
  }
  
}

void UKF::Prediction(double delta_t) {

  // STEP 1: Generate sigma points //

  // Create augmented state mean vector
  VectorXd x_aug = VectorXd(n_aug_); 
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Create augmented state covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_); 
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  MatrixXd L = P_aug.llt().matrixL(); // sqrt of P (lower triangular matrix)

  // Create augmented sigma points matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_aug.fill(0.0);
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i){
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  // debug print
  if(debug_flg) std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  // STEP 2: Predict sigma points //

  for (int i=0; i<2*n_aug_+1; ++i){
    // extract augmented states
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state process model
    double px_p, py_p;
    if (fabs(yawd) > 0.00001){ // avoid division by zero
      px_p = p_x + v/yawd * (sin(yaw+yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw+yawd*delta_t));
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);      
    }
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add process noise contribution
    px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_p = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // assign predicted sigma points
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // debug print
  if(debug_flg) std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;
  
  // STEP 3: Predict state mean and covariance //
  
  // Predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i){
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // Predicted state covariance matrix
  P_.fill(0);
  for (int i = 0; i < 2*n_aug_+1; ++i){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI){
      x_diff(3) -= 2.0*M_PI;
    }
    while (x_diff(3) < -M_PI){
      x_diff(3) += 2.0*M_PI;
    }
    P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
  }

  // debug print
  if(debug_flg) std::cout << "Predicted state = " << std::endl << x_ << std::endl;
  if(debug_flg) std::cout << "Predicted covariance = " << std::endl << P_ << std::endl;

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * Use lidar data to update the object's estimated stae, covariance, and the NIS.
   */

  // STEP 1: Measurement Prediction

  // measurment dimension (px, py)
  int n_z = 2;

  // sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // apply measurement model on sigma points
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    // innovation covariance
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
  S = S + R;

  // debug print
  if(debug_flg) std::cout << "Predicted measurement: " << std::endl << z_pred << std::endl;
  if(debug_flg) std::cout << "Innovation covariance: " << std::endl << P_ << std::endl;

  // STEP 2: State Update

  // Unpack measurement
  VectorXd z  = VectorXd(n_z);
  z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);

  // Cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3) > M_PI){
      x_diff(3) -= 2.0*M_PI;
    }
    while (x_diff(3) < -M_PI){
      x_diff(3) += 2.0*M_PI;
    }
    
    // cross correlation
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain
  MatrixXd K = Tc * S.inverse();
  if(debug_flg) std::cout << "Kalman gain: " << std::endl << K << std::endl;

  // residual
  VectorXd z_diff = z - z_pred;
  if(debug_flg) std::cout << "Measurement residual: " << std::endl << z_diff << std::endl;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // Compute normalized innovation squared (NIS)
  nis_lidar_ = z_diff.transpose()*S.inverse()*z_diff;

  // print result
  if(debug_flg) std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  if(debug_flg) std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  if(debug_flg) std::cout << "Lidar NIS: " << std::endl << nis_lidar_ << std::endl;
  if(print_flg){
		FILE *fpt;
		fpt = fopen("nis_out.csv", "a");
	  fprintf(fpt, "%d, %f, %f\n", 0, meas_package.timestamp_/1e6, nis_lidar_);
		fclose(fpt);
	}

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * Use Radar data to update the object's estimated stae, covariance, and NIS.
   */
  
  // STEP 1: Measurement Prediction

  // measurment dimension (r, phi, r_dot)
  int n_z = 3;

  // sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // apply measurement model on sigma points
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                     // r
    Zsig(1,i) = atan2(p_y,p_x);                              // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y); // r_dot
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI){
      z_diff(1) -= 2.0*M_PI;
    }
    while (z_diff(1) < -M_PI){
      z_diff(1) += 2.0*M_PI;
    }
    
    // innovation covariance
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;
  S = S + R;

  // debug print
  if(debug_flg) std::cout << "Predicted measurement: " << std::endl << z_pred << std::endl;
  if(debug_flg) std::cout << "Innovation covariance: " << std::endl << P_ << std::endl;

  // STEP 2: State Update

  // Unpack measurement
  VectorXd z = meas_package.raw_measurements_;

  // Cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI){
     z_diff(1) -= 2.0*M_PI;
    } 
    while (z_diff(1) < -M_PI) {
      z_diff(1) += 2.0*M_PI;
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3) > M_PI){
      x_diff(3) -= 2.0*M_PI;
    }
    while (x_diff(3) < -M_PI){
      x_diff(3) += 2.0*M_PI;
    }
    
    // cross correlation
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1) > M_PI){
    z_diff(1) -= 2.0*M_PI;
  }
  while (z_diff(1) < -M_PI){
    z_diff(1) += 2.0*M_PI;
  }

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // Compute normalized innovation squared (NIS)
  nis_radar_ = z_diff.transpose()*S.inverse()*z_diff;

  // print result
  if(debug_flg) std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  if(debug_flg) std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  if(debug_flg) std::cout << "Radar NIS: " << std::endl << nis_radar_ << std::endl;
  if(print_flg){
		FILE *fpt;
		fpt = fopen("nis_out.csv", "a");
	  fprintf(fpt, "%d, %f, %f\n", 1, meas_package.timestamp_/1e6, nis_radar_);
		fclose(fpt);
	}
 
}