 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  �ݻ�������
    %  状�?方程：x(:,k+1) = F * x(:,k) +[sqrt(Q) * randn;0]; 
    %  观测方程：z(k+1) = atan(0.1 * x(1,k+1)) + sqrt(R) * randn; 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc; clear all;close all;
    n=501;
    tf = 500;                                     % 模拟长度 
    x=zeros(2,n);
    z=zeros(1,n);
    x(:,1) =[1;0.1];                              % 初始状�? 
    x_ckf=zeros(2,n);
    % x_estimate(:,1) = [1;0.1];                  %状�?的估�?
    x_ckf(:,1)=[1;0.1];
    % e_x_estimate = x_estimate(:,1);             %EKF的初始估�?
    xhat=x_ckf(:,1);
    x_e_error=zeros(1,n);
    x_c_error=zeros(1,n);
    z_e_error=zeros(1,n);
    z_c_error=zeros(1,n);
    Q = 0.0001;                                    % 过程状�?协方�?
    R = 0.16;                                      % 测量噪声协方�?
    P =[0.0099,0;0,0.0001];                        %初始估计方差
    Pplus=P;
    F=[1,1;0,1];
    Gamma=[0.5;1];
    w=0.25;  
    kesi=sqrt(2)*[1,0,-1,0;0,1,0,-1];
    for k = 1 : tf 
        % 模拟系统 
       x(:,k+1) = F * x(:,k) + Gamma * sqrt(Q) * randn;      %状�?�?
       %x(:,k+1) = F * x(:,k) +[sqrt(Q) * randn;0]; 
       z(k+1) = atan(0.1 * x(1,k+1)) + sqrt(R) * randn;      %观测�?
    end;
    for k = 1 : tf 
        %Cubature卡尔曼滤波器
        %%%%%�?）求协方差矩阵平方根
        S=chol(Pplus,'lower');
        %%%%%�?）计算求容积�?
        rjpoint(:,1)=S*kesi(:,1)+xhat;
        rjpoint(:,2)=S*kesi(:,2)+xhat;
        rjpoint(:,3)=S*kesi(:,3)+xhat;
        rjpoint(:,4)=S*kesi(:,4)+xhat;
        %%%%%�?）传播求容积�?
        Xminus(:,1)=F*rjpoint(:,1);                           %容积点经过非线�?函数后的�?
        Xminus(:,2)=F*rjpoint(:,2);
        Xminus(:,3)=F*rjpoint(:,3); 
        Xminus(:,4)=F*rjpoint(:,4); 
        %%%%�?）状态预�?
        xminus=w*Xminus(:,1)+w*Xminus(:,2)+w*Xminus(:,3)+w*Xminus(:,4);
        %%%%(5)状�?预测协方差阵
        Pminus=w*(Xminus(:,1)*Xminus(:,1)'+Xminus(:,2)*Xminus(:,2)'+Xminus(:,3)*Xminus(:,3)'+Xminus(:,4)*Xminus(:,4)')-xminus*xminus'+Gamma * Q* Gamma';
       %Pminus=w*(Xminus(:,1)*Xminus(:,1)'+Xminus(:,2)*Xminus(:,2)'+Xminus(:,3)*Xminus(:,3)'+Xminus(:,4)*Xminus(:,4)')-xminus*xminus'+[Q,0;0,0]; 
       %%%%观测更新
        %%%%%�?）矩阵分�?
        Sminus=chol(Pminus,'lower');
        %%%%%�?）计算求容积�?
        rjpoint1(:,1)=Sminus*kesi(:,1)+xminus;
        rjpoint1(:,2)=Sminus*kesi(:,2)+xminus;
        rjpoint1(:,3)=Sminus*kesi(:,3)+xminus;
        rjpoint1(:,4)=Sminus*kesi(:,4)+xminus;
        %%%%%�?）传播求容积�?
        Z(1)=atan(0.1*rjpoint1(1,1));
        Z(2)=atan(0.1*rjpoint1(1,2));
        Z(3)=atan(0.1*rjpoint1(1,3));
        Z(4)=atan(0.1*rjpoint1(1,4));
       % Z(:,4)=[atan(0.1*rjpoint1(1,4));0];
        %%%%%%%�?）观测预�?
        zhat=w*(Z(1)+Z(2)+Z(3)+Z(4));
        %%%%(5)观测预测协方差阵
        %Pzminus=w*(Z(:,1)*Z(:,1)'+Z(:,2)*Z(:,2)'+Z(:,3)*Z(:,3)'+Z(:,4)*Z(:,4)')-zhat*zhat'+[R,0;0,Q];
        Pzminus=w*(Z(1)^2+Z(2)^2+Z(3)^2+Z(4)^2)-zhat^2+R;
        %%%%(6)互协方差�?
        Pxzminus=w*(rjpoint1(:,1)*Z(1)+rjpoint1(:,2)*Z(2)+rjpoint1(:,3)*Z(3)+rjpoint1(:,4)*Z(4))-xminus*zhat;
        %%%%(7)计算卡尔曼增�?
        K=Pxzminus/Pzminus;
        %%%%(8)状�?更新
        xhat=xminus+K*(z(k+1)-zhat);
        %%%%(9)状�?协方差矩阵更�?
        Pplus=Pminus-K*Pzminus*K';
        
        x_ckf(:,k+1)=xhat;
    end
        t = 0 : tf;
    figure;
    plot(t,x(1,:),'k.',t,x_ckf(1,:),'g');
    legend('真实�?,'CKF估计�?);
