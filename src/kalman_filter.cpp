#include "kalman_filter.h"
#include <math.h>
KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

#define PI 3.14159265

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    //Convert from cartezian x_ to polar radar measurement (non linear conversion) 
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
    float range = sqrt(px*px + py*py);
    VectorXd z_pred = VectorXd(3);
    if (abs(px) > 0.001) 
    {
        z_pred(0) = range;
        z_pred(1) = atan(py/px);
        z_pred(2) = (px*vx + py*vy)/range;
    } else {
        z_pred(0) = py;
        z_pred(1) = PI/2;
        z_pred(2) = vy;
    }
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
    
}

using namespace std;
#include <iostream>

//Calculates the Normalized Innovation Squared vector
void KalmanFilter::CalculateNIS(const VectorXd& gt) {
  VectorXd RES = H_*(x_ - gt);
  MatrixXd RESt = RES.transpose();
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  NIS_ = RESt * Si * RES;
}

