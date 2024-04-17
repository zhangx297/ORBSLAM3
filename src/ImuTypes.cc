/**
* This file is part of ORB-SLAM3
*
* Copyright (C) 2017-2021 Carlos Campos, Richard Elvira, Juan J. Gómez Rodríguez, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
* Copyright (C) 2014-2016 Raúl Mur-Artal, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
*
* ORB-SLAM3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with ORB-SLAM3.
* If not, see <http://www.gnu.org/licenses/>.
*/

#include "ImuTypes.h"
#include "Converter.h"

#include "GeometricTools.h"

#include<iostream>
#include <fstream>

namespace ORB_SLAM3
{

namespace IMU
{

const float eps = 1e-4;

// 旋转矩阵的奇异值分解
Eigen::Matrix3f NormalizeRotation(const Eigen::Matrix3f &R){
    Eigen::JacobiSVD<Eigen::Matrix3f> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);  //奇异值分解的雅可比方法
    return svd.matrixU() * svd.matrixV().transpose();
}

// 右乘雅可比
Eigen::Matrix3f RightJacobianSO3(const float &x, const float &y, const float &z)
{
    Eigen::Matrix3f I;
    I.setIdentity();  //单位阵
    const float d2 = x*x+y*y+z*z;
    const float d = sqrt(d2);
    Eigen::Vector3f v;
    v << x, y, z;
    Eigen::Matrix3f W = Sophus::SO3f::hat(v); //3x1的向量转换为3x3的反对称矩阵
    if(d<eps) {
        return I;
    }
    else {
        return I - W*(1.0f-cos(d))/d2 + W*W*(d-sin(d))/(d2*d); //雅可比
    }
}

Eigen::Matrix3f RightJacobianSO3(const Eigen::Vector3f &v)
{
    return RightJacobianSO3(v(0),v(1),v(2));
}

// 右乘雅可比求逆
Eigen::Matrix3f InverseRightJacobianSO3(const float &x, const float &y, const float &z)
{
    Eigen::Matrix3f I;
    I.setIdentity();
    const float d2 = x*x+y*y+z*z;
    const float d = sqrt(d2);
    Eigen::Vector3f v;
    v << x, y, z;
    Eigen::Matrix3f W = Sophus::SO3f::hat(v);

    if(d<eps) {
        return I;
    }
    else {
        return I + W/2 + W*W*(1.0f/d2 - (1.0f+cos(d))/(2.0f*d*sin(d)));
    }
}

Eigen::Matrix3f InverseRightJacobianSO3(const Eigen::Vector3f &v)
{
    return InverseRightJacobianSO3(v(0),v(1),v(2));
}

// 角速度计算deltaR，的右乘雅克比rightJ
IntegratedRotation::IntegratedRotation(const Eigen::Vector3f &angVel, const Bias &imuBias, const float &time) {
    //角速度减去偏置，再乘上时间，就是该时间内的旋转，xyz构成旋转向量   
    const float x = (angVel(0)-imuBias.bwx)*time;
    const float y = (angVel(1)-imuBias.bwy)*time;
    const float z = (angVel(2)-imuBias.bwz)*time;
    //计算旋转向量的模值
    const float d2 = x*x+y*y+z*z;
    const float d = sqrt(d2);

    Eigen::Vector3f v;
    v << x, y, z;
    //旋转向量（so3）写成反对程形式
    Eigen::Matrix3f W = Sophus::SO3f::hat(v);
    if(d<eps)  ////旋转比较小，旋转向量到旋转矩阵的指数映射采用一阶近似
    {   
        //forster 经典预积分公式（4）
        deltaR = Eigen::Matrix3f::Identity() + W;
        //小量时，右扰动  Jr = I
        rightJ = Eigen::Matrix3f::Identity();
    }
    else
    {
        //forster 经典预积分公式（3）// deltaR
        deltaR = Eigen::Matrix3f::Identity() + W*sin(d)/d + W*W*(1.0f-cos(d))/d2;
        //forster 经典预积分公式（8）// 右乘雅可比
        rightJ = Eigen::Matrix3f::Identity() - W*(1.0f-cos(d))/d2 + W*W*(d-sin(d))/(d2*d);
    }
}

// 构造预积分器
Preintegrated::Preintegrated(const Bias &b_, const Calib &calib)
{
    Nga = calib.Cov; //imu数据的高斯噪声协方差(6*6)
    NgaWalk = calib.CovWalk; //imu的bias随机游走的协方差(6*6)
    Initialize(b_);
}

// Copy constructor
Preintegrated::Preintegrated(Preintegrated* pImuPre): dT(pImuPre->dT),C(pImuPre->C), Info(pImuPre->Info),
     Nga(pImuPre->Nga), NgaWalk(pImuPre->NgaWalk), b(pImuPre->b), dR(pImuPre->dR), dV(pImuPre->dV),
    dP(pImuPre->dP), JRg(pImuPre->JRg), JVg(pImuPre->JVg), JVa(pImuPre->JVa), JPg(pImuPre->JPg), JPa(pImuPre->JPa),
    avgA(pImuPre->avgA), avgW(pImuPre->avgW), bu(pImuPre->bu), db(pImuPre->db), mvMeasurements(pImuPre->mvMeasurements)
{

}

void Preintegrated::CopyFrom(Preintegrated* pImuPre)
{
    dT = pImuPre->dT;
    C = pImuPre->C;
    Info = pImuPre->Info;
    Nga = pImuPre->Nga;
    NgaWalk = pImuPre->NgaWalk;
    b.CopyFrom(pImuPre->b);
    dR = pImuPre->dR;
    dV = pImuPre->dV;
    dP = pImuPre->dP;
    JRg = pImuPre->JRg;
    JVg = pImuPre->JVg;
    JVa = pImuPre->JVa;
    JPg = pImuPre->JPg;
    JPa = pImuPre->JPa;
    avgA = pImuPre->avgA;
    avgW = pImuPre->avgW;
    bu.CopyFrom(pImuPre->bu);
    db = pImuPre->db;
    mvMeasurements = pImuPre->mvMeasurements;
}

// 预积分参数定义
void Preintegrated::Initialize(const Bias &b_)
{
    dR.setIdentity();  //旋转预积分初始值为单位阵
    dV.setZero(); //速度预积分初始为0
    dP.setZero(); //位置预积分初始为0
    JRg.setZero(); //旋转预积分对delta(bg)的雅克比，在bias改变后预积分量一阶近似更新模型中使用
    JVg.setZero(); //速度预积分对delta(bg)的雅克比
    JVa.setZero(); //速度预积分对delta(ba)的雅克比
    JPg.setZero(); //位置预积分对delta(bg)的雅克比
    JPa.setZero(); //位置预积分对delta(ba)的雅克比
    C.setZero(); //协方差传递所需的A,B矩阵，A为左上角9*9，B为右下角6*6
    Info.setZero(); //信息矩阵(协方差矩阵的逆)
    db.setZero(); //bias的变化量
    b=b_; //起初偏差
    bu=b_; //updated bias
    avgA.setZero(); //平均加速度，用于判断加速度是否变化
    avgW.setZero(); //平均角速度
    dT=0.0f; //时间间隔
    mvMeasurements.clear();
}

// 重建预积分
void Preintegrated::Reintegrate()
{
    std::unique_lock<std::mutex> lock(mMutex);
    const std::vector<integrable> aux = mvMeasurements; //新测量值？
    Initialize(bu); //bu,更新的偏差
    for(size_t i=0;i<aux.size();i++)
        IntegrateNewMeasurement(aux[i].a,aux[i].w,aux[i].t);
}

// 3. 利用新的imu数据更新预积分量
void Preintegrated::IntegrateNewMeasurement(const Eigen::Vector3f &acceleration, const Eigen::Vector3f &angVel, const float &dt)
{
    //将imu数据构造成一个integrable结构体，保存到mvMeasurements中
    mvMeasurements.push_back(integrable(acceleration,angVel,dt));

    // Position is updated firstly, as it depends on previously computed velocity and rotation.
    // Velocity is updated secondly, as it depends on previously computed rotation.
    // Rotation is the last to be updated.
    // Matrices to compute covariance
    //计算协方差传递所需的A、B矩阵，下面的计算详见forster论文附录A.7~A.9
    Eigen::Matrix<float,9,9> A;
    A.setIdentity();
    Eigen::Matrix<float,9,6> B;
    B.setZero();

    // 减去偏置后的加速度、角速度
    Eigen::Vector3f acc, accW;
    acc << acceleration(0)-b.bax, acceleration(1)-b.bay, acceleration(2)-b.baz;
    accW << angVel(0)-b.bwx, angVel(1)-b.bwy, angVel(2)-b.bwz;

    // 计算平均加速度和角速度
    avgA = (dT*avgA + dR*acc*dt)/(dT+dt);
    avgW = (dT*avgW + accW*dt)/(dT+dt);

    //更新P、V的预积分量，forster论文公式(26)
    // Update delta position dP and velocity dV (rely on no-updated delta rotation)
    dP = dP + dV*dt + 0.5f*dR*acc*dt*dt;  
    dV = dV + dR*acc*dt;

    // Compute velocity and position parts of matrices A and B (rely on non-updated delta rotation)
    // 根据η_ij = A * η_i,j-1 + B_j-1 * η_j-1中的Ａ矩阵和Ｂ矩阵对速度和位移进行更新 （62）
    Eigen::Matrix<float,3,3> Wacc = Sophus::SO3f::hat(acc); //acc加速度反对称矩阵

    A.block<3,3>(3,0) = -dR*dt*Wacc;
    A.block<3,3>(6,0) = -0.5f*dR*dt*dt*Wacc;
    A.block<3,3>(6,3) = Eigen::DiagonalMatrix<float,3>(dt, dt, dt);
    B.block<3,3>(3,3) = dR*dt;
    B.block<3,3>(6,3) = 0.5f*dR*dt*dt;

    // 更新后的dP, dV更新对bias的jacobian
    // 递推公式的推导与上文中dP, dV的更新类似，都是将整个求和项分为i~j-2与j-1两部分
    // Update position and velocity jacobians wrt bias correction
    JPa = JPa + JVa*dt -0.5f*dR*dt*dt;
    JPg = JPg + JVg*dt -0.5f*dR*dt*dt*Wacc*JRg;
    JVa = JVa - dR*dt;
    JVg = JVg - dR*dt*Wacc*JRg;

    //更新dR并更新A，B矩阵中涉及到更新后的dR的部分
    // Update delta rotation
    //对角速度进行积分，得到旋转变化量
    IntegratedRotation dRi(angVel,b,dt);
    // 旧的旋转预积分量乘上旋转变化量，归一化使其符合旋转矩阵的格式
    dR = NormalizeRotation(dR*dRi.deltaR);

    // Compute rotation parts of matrices A and B
    // 补充AB矩阵
    A.block<3,3>(0,0) = dRi.deltaR.transpose();
    B.block<3,3>(0,0) = dRi.rightJ*dt;

    // Update covariance（63）
    C.block<9,9>(0,0) = A * C.block<9,9>(0,0) * A.transpose() + B*Nga*B.transpose();
    // 这一部分最开始是0矩阵，随着积分次数增加，每次都加上随机游走，偏置的信息矩阵
    C.block<6,6>(9,9) += NgaWalk;

    // Update rotation jacobian wrt bias correction
    // 更新dR对陀螺仪bias的jacobian ∂ΔRij/∂bg = (ΔRjj-1) * ∂ΔRij-1/∂bg - Jr(j-1)*t
    // 这里必须先更新dRi才能更新到这个值
    JRg = dRi.deltaR.transpose()*JRg - dRi.rightJ*dt;

    // Total integrated time
    dT += dt;
}

//融合两个预积分，发生在删除关键帧的时候，3帧变2帧，需要把两段预积分融合
void Preintegrated::MergePrevious(Preintegrated* pPrev)
{
    if (pPrev==this)
        return;

    std::unique_lock<std::mutex> lock1(mMutex);
    std::unique_lock<std::mutex> lock2(pPrev->mMutex);
    Bias bav;
    bav.bwx = bu.bwx;
    bav.bwy = bu.bwy;
    bav.bwz = bu.bwz;
    bav.bax = bu.bax;
    bav.bay = bu.bay;
    bav.baz = bu.baz;

    const std::vector<integrable > aux1 = pPrev->mvMeasurements;
    const std::vector<integrable> aux2 = mvMeasurements;

    Initialize(bav);
    for(size_t i=0;i<aux1.size();i++)
        IntegrateNewMeasurement(aux1[i].a,aux1[i].w,aux1[i].t);
    for(size_t i=0;i<aux2.size();i++)
        IntegrateNewMeasurement(aux2[i].a,aux2[i].w,aux2[i].t);

}

//设定新偏差
void Preintegrated::SetNewBias(const Bias &bu_)
{
    std::unique_lock<std::mutex> lock(mMutex);
    bu = bu_; //偏差

    db(0) = bu_.bwx-b.bwx;
    db(1) = bu_.bwy-b.bwy;
    db(2) = bu_.bwz-b.bwz;
    db(3) = bu_.bax-b.bax;
    db(4) = bu_.bay-b.bay;
    db(5) = bu_.baz-b.baz;
}

// 获得当前偏置与输入偏置的改变量
IMU::Bias Preintegrated::GetDeltaBias(const Bias &b_)
{
    std::unique_lock<std::mutex> lock(mMutex);
    return IMU::Bias(b_.bax-b.bax,b_.bay-b.bay,b_.baz-b.baz,b_.bwx-b.bwx,b_.bwy-b.bwy,b_.bwz-b.bwz);
}

// 根据新的偏置计算新的dR （44）
Eigen::Matrix3f Preintegrated::GetDeltaRotation(const Bias &b_)
{
    std::unique_lock<std::mutex> lock(mMutex);
    // 计算偏置的变化量
    Eigen::Vector3f dbg;
    dbg << b_.bwx-b.bwx,b_.bwy-b.bwy,b_.bwz-b.bwz;   //角速度的偏差
    // std::cout << NormalizeRotation(Sophus::SO3f::exp(JRg * dbg).matrix()) << std::endl;
    // std::cout << "--------------------------------------" << std::endl;

    
    // std::ofstream file("bR.txt", std::ofstream::app);
    // if(file.is_open()) {
    
    //     file << NormalizeRotation(Sophus::SO3f::exp(JRg * dbg).matrix()) << std::endl;
    //     file.close();
    // } else {
    //     std::cout << "文件打开失败." << std::endl;
    // }


    return NormalizeRotation(dR * Sophus::SO3f::exp(JRg * dbg).matrix());  //JRg旋转预积分对delta(bg)的雅克比，dgb偏置的变化
}

// 根据新的偏置计算新的dV (44)
Eigen::Vector3f Preintegrated::GetDeltaVelocity(const Bias &b_)
{
    std::unique_lock<std::mutex> lock(mMutex);
    Eigen::Vector3f dbg, dba;
    dbg << b_.bwx-b.bwx,b_.bwy-b.bwy,b_.bwz-b.bwz; //角速度的偏差
    dba << b_.bax-b.bax,b_.bay-b.bay,b_.baz-b.baz; //加速度的偏差
    Eigen::Vector3f a;
    a = JVg * dbg + JVa * dba;
    // std::cout << a[0] << " " << a[1] << " " << a[2] << std::endl; // JV* Eigen::Matrix3f   bax, bay, baz
    // std::cout << "--------------------------------------" << std::endl;


    // std::ofstream file("bV.txt", std::ofstream::app);
    // if(file.is_open()) {
    //     file << a[0] << " " << a[1] << " " << a[2] << std::endl;
    //     file.close();
    // } else {
    //     std::cout << "文件打开失败." << std::endl;
    // }


    return dV + JVg * dbg + JVa * dba;  //一阶近似更新
}

// 根据新的偏置计算新的dP (44)
Eigen::Vector3f Preintegrated::GetDeltaPosition(const Bias &b_)
{
    std::unique_lock<std::mutex> lock(mMutex);
    Eigen::Vector3f dbg, dba;
    dbg << b_.bwx-b.bwx,b_.bwy-b.bwy,b_.bwz-b.bwz; //角速度的偏差
    dba << b_.bax-b.bax,b_.bay-b.bay,b_.baz-b.baz; //加速度的偏差
    Eigen::Vector3f a;
    a = JPg * dbg + JPa * dba;
    // std::cout << a[0] << " " << a[1] << " " << a[2] << std::endl;
    // std::cout << "--------------------------------------" << std::endl;


    // std::ofstream file("bP.txt", std::ofstream::app);
    // if(file.is_open()) {
    //     file << a[0] << " " << a[1] << " " << a[2] << std::endl;
    //     file.close();
    // } else {
    //     std::cout << "文件打开失败." << std::endl;
    // }


    return dP + JPg * dbg + JPa * dba;  //一阶近似更新
}

// 根据总偏置计算更新的dR
Eigen::Matrix3f Preintegrated::GetUpdatedDeltaRotation()
{
    std::unique_lock<std::mutex> lock(mMutex);
    return NormalizeRotation(dR * Sophus::SO3f::exp(JRg*db.head(3)).matrix());  //db偏差的变换量
}

// 根据总偏置计算更新的dV
Eigen::Vector3f Preintegrated::GetUpdatedDeltaVelocity()
{
    std::unique_lock<std::mutex> lock(mMutex);
    return dV + JVg * db.head(3) + JVa * db.tail(3);
}

// 根据总偏置计算更新的dP
Eigen::Vector3f Preintegrated::GetUpdatedDeltaPosition()
{
    std::unique_lock<std::mutex> lock(mMutex);
    return dP + JPg*db.head(3) + JPa*db.tail(3);
}

Eigen::Matrix3f Preintegrated::GetOriginalDeltaRotation() {
    std::unique_lock<std::mutex> lock(mMutex);
    return dR;
}

Eigen::Vector3f Preintegrated::GetOriginalDeltaVelocity() {
    std::unique_lock<std::mutex> lock(mMutex);
    return dV;
}

Eigen::Vector3f Preintegrated::GetOriginalDeltaPosition()
{
    std::unique_lock<std::mutex> lock(mMutex);
    return dP;
}

Bias Preintegrated::GetOriginalBias()
{
    std::unique_lock<std::mutex> lock(mMutex);
    return b;
}

Bias Preintegrated::GetUpdatedBias()
{
    std::unique_lock<std::mutex> lock(mMutex);
    return bu;
}

Eigen::Matrix<float,6,1> Preintegrated::GetDeltaBias()
{
    std::unique_lock<std::mutex> lock(mMutex);
    return db;
}

void Bias::CopyFrom(Bias &b)
{
    bax = b.bax;
    bay = b.bay;
    baz = b.baz;
    bwx = b.bwx;
    bwy = b.bwy;
    bwz = b.bwz;
}

std::ostream& operator<< (std::ostream &out, const Bias &b)
{
    if(b.bwx>0)
        out << " ";
    out << b.bwx << ",";
    if(b.bwy>0)
        out << " ";
    out << b.bwy << ",";
    if(b.bwz>0)
        out << " ";
    out << b.bwz << ",";
    if(b.bax>0)
        out << " ";
    out << b.bax << ",";
    if(b.bay>0)
        out << " ";
    out << b.bay << ",";
    if(b.baz>0)
        out << " ";
    out << b.baz;

    return out;
}

void Calib::Set(const Sophus::SE3<float> &sophTbc, const float &ng, const float &na, const float &ngw, const float &naw) {
    mbIsSet = true;
    const float ng2 = ng*ng;
    const float na2 = na*na;
    const float ngw2 = ngw*ngw;
    const float naw2 = naw*naw;

    // Sophus/Eigen
    //相机到imu变换矩阵
    mTbc = sophTbc;
    //imu到相机的变换矩阵
    mTcb = mTbc.inverse();
    Cov.diagonal() << ng2, ng2, ng2, na2, na2, na2;  //na、ng、naw、ngw都是标定的imu参数
    CovWalk.diagonal() << ngw2, ngw2, ngw2, naw2, naw2, naw2;
}

//IMU标定，相机IMU转换矩阵以及噪声
Calib::Calib(const Calib &calib)
{
    mbIsSet = calib.mbIsSet;
    // Sophus/Eigen parameters
    mTbc = calib.mTbc;
    mTcb = calib.mTcb;
    Cov = calib.Cov;
    CovWalk = calib.CovWalk;
}

} //namespace IMU

} //namespace ORB_SLAM2
