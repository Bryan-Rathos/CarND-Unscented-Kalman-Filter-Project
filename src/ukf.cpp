#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

float NormalizeAngle(float phi_in);
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() 
{
  // Set to false initailly. Set to true in the first call to ProcessMeasurement() 
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  // Inititialize timestamp, in microseconds
  time_us_ = 0.0;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Predicted sigma point matrix
  Xsig_pred_ =  MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Weights vector
  weights_ = VectorXd(2 * n_aug_ + 1);

  // the current NIS for radar
  NIS_radar_ = 0.0;

  // the current NIS for laser
  NIS_laser_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) 
{
  /*TODO:Complete this function! Make sure you switch between lidar and radar measurements.*/

  /************************
  * Initialization 
  ************************/
  if((meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ ) || 
     (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ ) )
  {
    if(!is_initialized_)
    {
      x_ << 1, 1, 1, 1, 1;

      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

      if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true)
      {
        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);
      }

      if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true)
      {
        float rho      = meas_package.raw_measurements_(0);
        float phi      = meas_package.raw_measurements_(1);
        float rho _dot = meas_package.raw_measurements_(3);

        x_(0) = rho * cos(phi);
        x_(1) = rho * sin(phi);
      }

      // Record current timestamp
      time_us_ = meas_package.timestamp_;

      // Done initializing, use sensor measurements from next time stamp
      is_initialized_ = true;        
      return;
    }
    
    /**********************
    * Prediction
    **********************/
    //compute the time elapsed between the current and previous measurements
    float dt = (measurement_pack.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
    time_us_ = measurement_pack.timestamp_;

    Prediction(dt);

    /**********************
    * Update
    **********************/
    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      UpdateLidar(meas_package);
    }
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      UpdateRadar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
  TODO:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /************************
  * Generate sigma points
  ************************/
  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  
  //calculate square root of P
  MatrixXd A = P.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x; i++)
  {
    Xsig.col(i+1)     = x + sqrt(lambda+n_x) * A.col(i);
    Xsig.col(i+1+n_x) = x - sqrt(lambda+n_x) * A.col(i);
  }

  /****************************
  * Sigma points augmentation
  ****************************/
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  //set lambda for augmented sigma points
  lambda_ = 3 - n_aug_;

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_*std_a_;
  P_aug(6, 6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  /****************************
  * Sigma points prediction
  ****************************/

  //predict sigma points
  for(int i = 0; i < 2 * n_aug + 1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v   = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    
    //predicted state values
    double px_p, py_p;
    
    //avoid division by zero
    if(fabs(yawd) > 0.001)
    {
      px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t) );    
    }
    else
    {
      px_p = p_x + v * cos(yaw) * delta_t;
      py_p = p_y + v * sin(yaw) * delta_t;
    }
    
    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;
    
    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  /***************************************************************
  * Predict Mean and Covariance from the predeicted sigma points
  ***************************************************************/
  //set weights
  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i < 2 * n_aug_ + 1; i++)
  { 
    double weight = 0.5/(n_aug+lambda_);
    weights_(i) = weight;
  }

  //predict state mean
  x_.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);   //update state 
  }
  
  //predict state covariance matrix
  P_.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x;     //state difference
    if(x_diff(3) > M_PI || x_diff(3) < - M_PI)
    {
      x_diff(3) = NormalizeAngle(x_diff(3));       //angle normalization
    }
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();  //update covariance matrix
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /*TODO:
    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.
    You'll also need to calculate the lidar NIS. */
  // Set LIDAR measurement dmenson
  int n_z = 2;

  // Retrieve measurement and store in vector z
  VectorXd z = meas_package.raw_measurements_;

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ +1);

  //Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // Build measurement model
  for(int i=0; i<2*n_aug+1; i++)
  {
    // Extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
  
    // Measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;                     
  }

  // Mean predicted measurement
  z_pred.fill(0.0);
  for(int i = 0; i < 2 * n_aug_+1; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Predicted measurement covariance matrix S
  S.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << pow(std_laspx_,2),0,
       0,pow(std_laspy_,2);
  S = S + R;

  /***** UKF LIDAR Update ******/

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++)
  {
    // Residual (predicted sigma pts. and predicted measurement mean pts)
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual (Actual measurement and predicted measurement)
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // Calculate NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /*TODO:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the radar NIS.*/

  // Set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // Retrieve measurement and store in vector z
  VectorXd z = meas_package.raw_measurements_;

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // Build measurement model
  for(int i=0; i < 2 * n_aug_ + 1; i++)
  {
    // Extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
  
    double v_x = v * cos(yaw);
    double v_y = v * sin(yaw);
    
    // Transform sigma points into measurement space
    // Measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);             //rho
    Zsig(1,i) = atan2(p_y,p_x);                      //phi

    if(Zsig(0,1) != 0)
      Zsig(2,i) = (p_x*v_x + p_y*v_y)/Zsig(0,i);     //rho_dot
    else
      Zsig(2,i) = 0;
  }

  // Calculate mean predicted measurement z_pred
  z_pred.fill(0.0);
  for(int i = 0; i < 2 * n_aug_+1; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Predicted measurement covariance matrix S
  S.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // Angle normalization
    if(z_diff(1) > M_PI || z_diff(1) < - M_PI)
    {
      z_diff(1) = NormalizeAngle(z_diff(1));
    }

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << pow(std_radr_,2),0,0,
       0,pow(std_radphi_,2),0,
       0,0,pow(std_radrd_,2);

  S = S + R;

  /***** UKF Radar Update ******/

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++)
  {
    // Residual (predicted sigma pts. and predicted measurement mean pts)
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    if(z_diff(1) > M_PI || z_diff(1) < - M_PI)
    {
      z_diff(1) = NormalizeAngle(z_diff(1));
    }

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    if(x_diff(3) > M_PI || x_diff(3) < - M_PI)
    {
      x_diff(3) = NormalizeAngle(x_diff(3));
    }

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();


  //residual (Actual measurement and predicted measurement)
  VectorXd z_diff = z - z_pred;

  //angle normalization
  if(z_diff(1) > M_PI || z_diff(1) < - M_PI)
  {
    z_diff(1) = NormalizeAngle(z_diff(1));
  }
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // Calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}


float NormalizeAngle(float phi_in)
{
  float norm_phi = atan2(sin(phi_in),cos(phi_in));
  return norm_phi;
}
