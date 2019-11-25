 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  �ݻ�Kalman�˲�
    %  ״̬���̣�x(:,k+1) = F * x(:,k) +[sqrt(Q) * randn;0]; 
    %  �۲ⷽ�̣�z(k+1) = atan(0.1 * x(1,k+1)) + sqrt(R) * randn; 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc; clear all;close all;
    n=501;
    tf = 500;                                     % ģ�ⳤ�� 
    x=zeros(2,n);
    z=zeros(1,n);
    x(:,1) =[1;0.1];                              % ��ʼ״̬ 
    x_ckf=zeros(2,n);
    % x_estimate(:,1) = [1;0.1];                  %״̬�Ĺ���
    x_ckf(:,1)=[1;0.1];
    % e_x_estimate = x_estimate(:,1);             %EKF�ĳ�ʼ����
    xhat=x_ckf(:,1);
    x_e_error=zeros(1,n);
    x_c_error=zeros(1,n);
    z_e_error=zeros(1,n);
    z_c_error=zeros(1,n);
    Q = 0.0001;                                    % ����״̬Э���� 
    R = 0.16;                                      % ��������Э���� 
    P =[0.0099,0;0,0.0001];                        %��ʼ���Ʒ���
    Pplus=P;
    F=[1,1;0,1];
    Gamma=[0.5;1];
    w=0.25;  
    kesi=sqrt(2)*[1,0,-1,0;0,1,0,-1];
    for k = 1 : tf 
        % ģ��ϵͳ 
       x(:,k+1) = F * x(:,k) + Gamma * sqrt(Q) * randn;      %״ֵ̬ 
       %x(:,k+1) = F * x(:,k) +[sqrt(Q) * randn;0]; 
       z(k+1) = atan(0.1 * x(1,k+1)) + sqrt(R) * randn;      %�۲�ֵ
    end;
    for k = 1 : tf 
        %Cubature�������˲���
        %%%%%��1����Э�������ƽ����
        S=chol(Pplus,'lower');
        %%%%%��2���������ݻ���
        rjpoint(:,1)=S*kesi(:,1)+xhat;
        rjpoint(:,2)=S*kesi(:,2)+xhat;
        rjpoint(:,3)=S*kesi(:,3)+xhat;
        rjpoint(:,4)=S*kesi(:,4)+xhat;
        %%%%%��3���������ݻ���
        Xminus(:,1)=F*rjpoint(:,1);                           %�ݻ��㾭�������Ժ������ֵ
        Xminus(:,2)=F*rjpoint(:,2);
        Xminus(:,3)=F*rjpoint(:,3); 
        Xminus(:,4)=F*rjpoint(:,4); 
        %%%%��4��״̬Ԥ��
        xminus=w*Xminus(:,1)+w*Xminus(:,2)+w*Xminus(:,3)+w*Xminus(:,4);
        %%%%(5)״̬Ԥ��Э������
        Pminus=w*(Xminus(:,1)*Xminus(:,1)'+Xminus(:,2)*Xminus(:,2)'+Xminus(:,3)*Xminus(:,3)'+Xminus(:,4)*Xminus(:,4)')-xminus*xminus'+Gamma * Q* Gamma';
       %Pminus=w*(Xminus(:,1)*Xminus(:,1)'+Xminus(:,2)*Xminus(:,2)'+Xminus(:,3)*Xminus(:,3)'+Xminus(:,4)*Xminus(:,4)')-xminus*xminus'+[Q,0;0,0]; 
       %%%%�۲����
        %%%%%��1������ֽ�
        Sminus=chol(Pminus,'lower');
        %%%%%��2���������ݻ���
        rjpoint1(:,1)=Sminus*kesi(:,1)+xminus;
        rjpoint1(:,2)=Sminus*kesi(:,2)+xminus;
        rjpoint1(:,3)=Sminus*kesi(:,3)+xminus;
        rjpoint1(:,4)=Sminus*kesi(:,4)+xminus;
        %%%%%��3���������ݻ���
        Z(1)=atan(0.1*rjpoint1(1,1));
        Z(2)=atan(0.1*rjpoint1(1,2));
        Z(3)=atan(0.1*rjpoint1(1,3));
        Z(4)=atan(0.1*rjpoint1(1,4));
       % Z(:,4)=[atan(0.1*rjpoint1(1,4));0];
        %%%%%%%��4���۲�Ԥ��
        zhat=w*(Z(1)+Z(2)+Z(3)+Z(4));
        %%%%(5)�۲�Ԥ��Э������
        %Pzminus=w*(Z(:,1)*Z(:,1)'+Z(:,2)*Z(:,2)'+Z(:,3)*Z(:,3)'+Z(:,4)*Z(:,4)')-zhat*zhat'+[R,0;0,Q];
        Pzminus=w*(Z(1)^2+Z(2)^2+Z(3)^2+Z(4)^2)-zhat^2+R;
        %%%%(6)��Э������
        Pxzminus=w*(rjpoint1(:,1)*Z(1)+rjpoint1(:,2)*Z(2)+rjpoint1(:,3)*Z(3)+rjpoint1(:,4)*Z(4))-xminus*zhat;
        %%%%(7)���㿨��������
        K=Pxzminus/Pzminus;
        %%%%(8)״̬����
        xhat=xminus+K*(z(k+1)-zhat);
        %%%%(9)״̬Э����������
        Pplus=Pminus-K*Pzminus*K';
        
        x_ckf(:,k+1)=xhat;
    end
        t = 0 : tf;
    figure;
    plot(t,x(1,:),'k.',t,x_ckf(1,:),'g');
    legend('��ʵֵ','CKF����ֵ');
