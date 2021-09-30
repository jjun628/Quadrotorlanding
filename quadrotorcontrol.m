%% Quadrotor UAV tracking simulation model / period: 21. 9.12. ~ / recent review: 21. 9.30.
%% Task
% 1. Add rotation estimation (Function)
% 2. Add disturbances (Wind)
% 3. Linearization, LQR controller

clear all; close all; clc

graphs = 1;
simulation = 1;


%% simulation time vector
tf = 10;
dt = 0.01;
t = 0:dt:tf;

%% variables
g = 9.81;  % gravitational accleration: m/s^2
m = 0.03;     % drone mass: kg
disturbance = [5; 0; 0];    % m/s,  10m/s = 36km/h
L = 0.046;    % arm length / meter
LL = [0, -L,  0,  L;      % drone arm vector in the body frame
      L,  0, -L,  0;
      0,  0,  0,  0];

eer = 0.05;        % moving platform edge size: 10cm
eeh = 0.005;       % moving platform height
bb= [ -eer,  eer,  eer, -eer, -eer,  eer,  eer, -eer;    % platform vector
      -eer, -eer,  eer,  eer, -eer, -eer,  eer,  eer;
       eeh,  eeh,  eeh,  eeh, -eeh, -eeh, -eeh, -eeh];  
   
%% mark vectors on the platform
mml= 0.04; 
mm = [  mml,    0, -mml,    0;
          0,  mml,    0, -mml;
        eeh,  eeh,  eeh,  eeh];

%% control gains
kp  = 100; kd  =  20;      % position controller
kp4 = 500; kp5 = 500; kp6 = 500; % attitude controller
kd4 = 10; kd5 = 10; kd6 = 10;

KF = 6.11e-8;    % aerodynamic force constant
KM = 1.5e-9;    % aerodynamic moment constant
gamma = KM/KF;

%% moment of inertia
% Ixx = 1/4*m*0.1^2 + 1/12*m*0.1^2; % drone body / meters
% Iyy = 1/4*m*0.1^2 + 1/12*m*0.1^2;
% Izz = 1/2*m*0.1^2;
Ixx = 1.43e-5; Iyy = 1.43e-5; Izz = 2.89e-5;
Ig = diag([Ixx,Iyy,Izz]);
Jr = 1/12*0.5*(4*L^2+0.05^2);   % rotor inertia

%% pre-allocate variables
x = zeros(6,size(t,2));  % state vector
xdot = zeros(6,size(t,2));  % derivative of state vector

%% platform state
x_dp = readmatrix('motion_platform.txt');  % dt=0.01, 10 seconds periodic motions
% vel_dp = [zeros(6,1) diff(x_dp,1,2)]/dt;
% acc_dp = [zeros(6,1) diff(vel_dp,1,2)]/dt;

%% initial & final states
altitude = 1;
edge = 0.2159;
x_i = [edge; edge; edge+altitude; 0; 0; 0];  % initial state
x_f = [edge; edge; edge+altitude; 0; 0; 0];  % final state
xdot_i = [0; 0; 0; 0; 0; 0];  % initial velocity
xdot_f = [0; 0; 0; 0; 0; 0];  % final velocity

%% rotation estimation values
cof = 0.5;  % coefficient factor   0.5 /m

%% Trajectory maker
k = t/tf;
x_d = x_i + (x_f-x_i) * (10*k.^3 - 15*k.^4 + 6*k.^5); % 5th order trajectory
% x_d=x_dp;
% vel_d = [zeros(6,1) diff(x_d,1,2)]/dt;
% acc_d = [zeros(6,1) diff(vel_d,1,2)]/dt;
% x_d=round(x_d,4);
% writematrix(x_d,'check.txt'); 

x(:,1) = x_i;
xdot(:,1) = xdot_i;
xddot = [zeros(6,1) diff(xdot,1,2)]/dt;
U = [];   % control vector

for i=1:length(t)
    %% Rotation matrix, XYZ Tait-Bryan
    phi=x(4,i); the=x(5,i); psi=x(6,i);
    c1=cos(phi); c2=cos(the); c3=cos(psi); s1=sin(phi); s2=sin(the); s3=sin(psi);
    R = [c2*c3,          -c2*s3,              s2;
         c1*s3+c3*s1*s2,  c1*c3-s1*s2*s3, -c2*s1;
         s1*s3-c1*c3*s2,  c3*s1+c1*s2*s3,  c1*c2]; 
    for k=1:4
        LLR(:,k,i) = R*LL(:,k);
    end
    
    %% platform rotation
    phip=x_dp(4,i); thep=x_dp(5,i); psip=x_dp(6,i);
    c1p=cos(phip); c2p=cos(thep); c3p=cos(psip); s1p=sin(phip); s2p=sin(thep); s3p=sin(psip);
    Rp = [c2p*c3p,                    -c2p*s3p,     s2p;
          c1p*s3p+c3p*s1p*s2p,  c1p*c3p-s1p*s2p*s3p, -c2p*s1p;
          s1p*s3p-c1p*c3p*s2p,  c3p*s1p+c1p*s2p*s3p,  c1p*c2p]; 
    for m=1:8, bbR(:,m,i) = Rp*bb(:,m); end   % platform edges
    for m=1:4, mmR(:,m,i) = Rp*mm(:,m); end   % mark vectors
    
    %% getting desired position from the camera(x, y values) w/ sensor noise
    noise = -0.0015 + (0.0015+0.0015)*rand(2,1);  % interval (a,b) with the formula a + (b-a).*rand(N,1)
    x_d(1:2,i) = x_dp(1:2,i) + noise;  % difference between quadrotor and platform center(x,y) except for z component
    vel_d = [zeros(6,1) diff(x_d,1,2)]/dt;
    acc_d = [zeros(6,1) diff(vel_d,1,2)]/dt;
    
    %% Camera image vector
    for m=1:4, mmC(:,m,i)=R*mmR(:,m,i)*cof*(x(3,i)-edge); end   % considering quadrotor rotation, coefficient factor
    mmC(3,:) = 0;

    %% desired roll and pitch
    x_d(4,i) = 1/g * (acc_d(1,i)*sin(x(6,i)) - acc_d(2,i)*cos(x(6,i)));
    x_d(5,i) = 1/g * (acc_d(1,i)*cos(x(6,i)) + acc_d(2,i)*sin(x(6,i)));

    errp = x_d(:,i)-x(:,i);
    errd = vel_d(:,i)-xdot(:,i);
    %% attitude controller
    U2 = Ig * [kp4*(errp(4)) + kd4*(errd(4));
               kp5*(errp(5)) + kd5*(errd(5));
               kp6*(errp(6)) + kd6*(errd(6))];
    
    %% position controller
    U1 = (m*(acc_d(3,i) - kp*errp(3) - kd*errd(3)) + m*g);
    
    %% Control input
    U=[U1;U2];
    U_rec(:,i)=U;
    
    %% Solving Dynamic equation
    [x(:,i+1), xdot(:,i+1)] = rk4_ode([x(:,i); xdot(:,i)],m,R,U,g,Ig,dt);
    xddot = [zeros(6,1) diff(xdot,1,2)]/dt;
end

%% quadrotor graphs
if graphs,
    figure;
    subplot(2,1,1); hold on; grid on
    plot(t,x(1,1:end-1),'r');plot(t, x_dp(1,1:end),'r--')
    plot(t,x(2,1:end-1),'b');plot(t, x_dp(2,1:end),'b--')
    plot(t,x(3,1:end-1),'g');plot(t, x_dp(3,1:end),'g--')
    subplot(2,1,2); hold on; grid on
    plot(t,x(4,1:end-1)*180/pi,'r');plot(t, x_d(4,1:end)*180/pi,'r--')
    plot(t,x(5,1:end-1)*180/pi,'b');plot(t, x_d(5,1:end)*180/pi,'b--')
    plot(t,x(6,1:end-1)*180/pi,'g');plot(t, x_d(6,1:end)*180/pi,'g--')
end

%% plotting
if simulation,
    figure; view(-70,30); grid on; hold on; j=1; xlabel('x');ylabel('y');zlabel('z'); xlim([edge-0.2,edge+0.2]),ylim([edge-0.2,edge+0.2]),zlim([0-0.1,edge+1.2])
    %% quadrotor part
    % plot3(x_d(1,:),x_d(2,:),x_d(3,:),'r');
    p1=plot3(x(1,j),x(2,j),x(3,j),'g*');
    p2=plot3([x(1,j)+LLR(1,1,j),x(1,j)+LLR(1,3,j)],[x(2,j)+LLR(2,1,j),x(2,j)+LLR(2,3,j)],[x(3,j)+LLR(3,1,j),x(3,j)+LLR(3,3,j)],'k');
    p3=plot3([x(1,j)+LLR(1,2,j),x(1,j)+LLR(1,4,j)],[x(2,j)+LLR(2,2,j),x(2,j)+LLR(2,4,j)],[x(3,j)+LLR(3,2,j),x(3,j)+LLR(3,4,j)],'k');
    %% platform part
%     pp1=plot3(x_dp(1,j),x_dp(2,j),x_dp(3,j),'b*','MarkerSize',1);
    % upper shpae
    platform(1) = plot3([x_dp(1,j)+bbR(1,1,j),x_dp(1,j)+bbR(1,2,j)],[x_dp(2,j)+bbR(2,1,j),x_dp(2,j)+bbR(2,2,j)],[x_dp(3,j)+bbR(3,1,j),x_dp(3,j)+bbR(3,2,j)],'b','LineWidth',1.2);
    platform(2) = plot3([x_dp(1,j)+bbR(1,3,j),x_dp(1,j)+bbR(1,2,j)],[x_dp(2,j)+bbR(2,3,j),x_dp(2,j)+bbR(2,2,j)],[x_dp(3,j)+bbR(3,3,j),x_dp(3,j)+bbR(3,2,j)],'b','LineWidth',1.2);
    platform(3) = plot3([x_dp(1,j)+bbR(1,3,j),x_dp(1,j)+bbR(1,4,j)],[x_dp(2,j)+bbR(2,3,j),x_dp(2,j)+bbR(2,4,j)],[x_dp(3,j)+bbR(3,3,j),x_dp(3,j)+bbR(3,4,j)],'b','LineWidth',1.2);
    platform(4) = plot3([x_dp(1,j)+bbR(1,1,j),x_dp(1,j)+bbR(1,4,j)],[x_dp(2,j)+bbR(2,1,j),x_dp(2,j)+bbR(2,4,j)],[x_dp(3,j)+bbR(3,1,j),x_dp(3,j)+bbR(3,4,j)],'b','LineWidth',1.2);
    % side shape
    platform(5) = plot3([x_dp(1,j)+bbR(1,1,j),x_dp(1,j)+bbR(1,5,j)],[x_dp(2,j)+bbR(2,1,j),x_dp(2,j)+bbR(2,5,j)],[x_dp(3,j)+bbR(3,1,j),x_dp(3,j)+bbR(3,5,j)],'b','LineWidth',1.2);
    platform(6) = plot3([x_dp(1,j)+bbR(1,6,j),x_dp(1,j)+bbR(1,2,j)],[x_dp(2,j)+bbR(2,6,j),x_dp(2,j)+bbR(2,2,j)],[x_dp(3,j)+bbR(3,6,j),x_dp(3,j)+bbR(3,2,j)],'b','LineWidth',1.2);
    platform(7) = plot3([x_dp(1,j)+bbR(1,3,j),x_dp(1,j)+bbR(1,7,j)],[x_dp(2,j)+bbR(2,3,j),x_dp(2,j)+bbR(2,7,j)],[x_dp(3,j)+bbR(3,3,j),x_dp(3,j)+bbR(3,7,j)],'b','LineWidth',1.2);
    platform(8) = plot3([x_dp(1,j)+bbR(1,8,j),x_dp(1,j)+bbR(1,4,j)],[x_dp(2,j)+bbR(2,8,j),x_dp(2,j)+bbR(2,4,j)],[x_dp(3,j)+bbR(3,8,j),x_dp(3,j)+bbR(3,4,j)],'b','LineWidth',1.2);
    % bottom shape
    platform(9) = plot3([x_dp(1,j)+bbR(1,6,j),x_dp(1,j)+bbR(1,5,j)],[x_dp(2,j)+bbR(2,6,j),x_dp(2,j)+bbR(2,5,j)],[x_dp(3,j)+bbR(3,6,j),x_dp(3,j)+bbR(3,5,j)],'b','LineWidth',1.2);
    platform(10)= plot3([x_dp(1,j)+bbR(1,6,j),x_dp(1,j)+bbR(1,7,j)],[x_dp(2,j)+bbR(2,6,j),x_dp(2,j)+bbR(2,7,j)],[x_dp(3,j)+bbR(3,6,j),x_dp(3,j)+bbR(3,7,j)],'b','LineWidth',1.2);
    platform(11)= plot3([x_dp(1,j)+bbR(1,8,j),x_dp(1,j)+bbR(1,7,j)],[x_dp(2,j)+bbR(2,8,j),x_dp(2,j)+bbR(2,7,j)],[x_dp(3,j)+bbR(3,8,j),x_dp(3,j)+bbR(3,7,j)],'b','LineWidth',1.2);
    platform(12)= plot3([x_dp(1,j)+bbR(1,8,j),x_dp(1,j)+bbR(1,5,j)],[x_dp(2,j)+bbR(2,8,j),x_dp(2,j)+bbR(2,5,j)],[x_dp(3,j)+bbR(3,8,j),x_dp(3,j)+bbR(3,5,j)],'b','LineWidth',1.2);
    %% mark part
    mark(1) = plot3([x_dp(1,j)+mmR(1,1,j),x_dp(1,j)+mmR(1,3,j)],[x_dp(2,j)+mmR(2,1,j),x_dp(2,j)+mmR(2,3,j)],[x_dp(3,j)+mmR(3,1,j),x_dp(3,j)+mmR(3,3,j)],'b','MarkerSize',8);
    mark(2) = plot3([x_dp(1,j)+mmR(1,2,j),x_dp(1,j)+mmR(1,4,j)],[x_dp(2,j)+mmR(2,2,j),x_dp(2,j)+mmR(2,4,j)],[x_dp(3,j)+mmR(3,2,j),x_dp(3,j)+mmR(3,4,j)],'b','MarkerSize',8); 
    mark(3) = plot3(x_dp(1,j)+mmR(1,1,j),x_dp(2,j)+mmR(2,1,j),x_dp(3,j)+mmR(3,1,j),'bo','MarkerSize',3); 
    mark(4) = plot3(x_dp(1,j)+mmR(1,2,j),x_dp(2,j)+mmR(2,2,j),x_dp(3,j)+mmR(3,2,j),'bo','MarkerSize',3); 
    mark(5) = plot3(x_dp(1,j)+mmR(1,3,j),x_dp(2,j)+mmR(2,3,j),x_dp(3,j)+mmR(3,3,j),'bo','MarkerSize',3); 
    mark(6) = plot3(x_dp(1,j)+mmR(1,4,j),x_dp(2,j)+mmR(2,4,j),x_dp(3,j)+mmR(3,4,j),'bo','MarkerSize',3); 
    mark(7) = plot3((x_dp(1,j)+mmR(1,1,j)+x_dp(1,j)+mmR(1,3,j))/2,(x_dp(2,j)+mmR(2,1,j)+x_dp(2,j)+mmR(2,3,j))/2,(x_dp(3,j)+mmR(3,1,j)+x_dp(3,j)+mmR(3,3,j))/2,'bo','MarkerSize',3); 
    %% Camera image part
    camera(1) = plot3([x_dp(1,j)+mmC(1,1,j),x_dp(1,j)+mmC(1,3,j)],[x_dp(2,j)+mmC(2,1,j),x_dp(2,j)+mmC(2,3,j)],[mmC(3,1,j)+x(3,j)-0.1,mmC(3,3,j)+x(3,j)-0.1],'r','MarkerSize',8);
    camera(2) = plot3([x_dp(1,j)+mmC(1,2,j),x_dp(1,j)+mmC(1,4,j)],[x_dp(2,j)+mmC(2,2,j),x_dp(2,j)+mmC(2,4,j)],[mmC(3,2,j)+x(3,j)-0.1,mmC(3,4,j)+x(3,j)-0.1],'r','MarkerSize',8); 
    camera(3) = plot3(x_dp(1,j)+mmC(1,1,j),x_dp(2,j)+mmC(2,1,j),mmC(3,1,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
    camera(4) = plot3(x_dp(1,j)+mmC(1,2,j),x_dp(2,j)+mmC(2,2,j),mmC(3,2,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
    camera(5) = plot3(x_dp(1,j)+mmC(1,3,j),x_dp(2,j)+mmC(2,3,j),mmC(3,3,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
    camera(6) = plot3(x_dp(1,j)+mmC(1,4,j),x_dp(2,j)+mmC(2,4,j),mmC(3,4,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
    camera(7) = plot3((x_dp(1,j)+mmC(1,1,j)+x_dp(1,j)+mmC(1,3,j))/2,(x_dp(2,j)+mmC(2,1,j)+x_dp(2,j)+mmC(2,3,j))/2,(mmC(3,1,j)+mmC(3,3,j)+x(3,j)-0.1+x(3,j)-0.1)/2,'ro','MarkerSize',3); 
    pause; delete(p1);delete(p2);delete(p3); delete(platform); delete(mark); delete(camera)
    for j=1:length(t)
        %% quadrotor part
        p1=plot3(x(1,j),x(2,j),x(3,j),'g*');
        p2=plot3([x(1,j)+LLR(1,1,j),x(1,j)+LLR(1,3,j)],[x(2,j)+LLR(2,1,j),x(2,j)+LLR(2,3,j)],[x(3,j)+LLR(3,1,j),x(3,j)+LLR(3,3,j)],'k');
        p3=plot3([x(1,j)+LLR(1,2,j),x(1,j)+LLR(1,4,j)],[x(2,j)+LLR(2,2,j),x(2,j)+LLR(2,4,j)],[x(3,j)+LLR(3,2,j),x(3,j)+LLR(3,4,j)],'k');
        %% platform part
%         pp1=plot3(x_dp(1,j),x_dp(2,j),x_dp(3,j),'b*','MarkerSize',1);
        % upper shpae
        platform(1) = plot3([x_dp(1,j)+bbR(1,1,j),x_dp(1,j)+bbR(1,2,j)],[x_dp(2,j)+bbR(2,1,j),x_dp(2,j)+bbR(2,2,j)],[x_dp(3,j)+bbR(3,1,j),x_dp(3,j)+bbR(3,2,j)],'b','LineWidth',1.2);
        platform(2) = plot3([x_dp(1,j)+bbR(1,3,j),x_dp(1,j)+bbR(1,2,j)],[x_dp(2,j)+bbR(2,3,j),x_dp(2,j)+bbR(2,2,j)],[x_dp(3,j)+bbR(3,3,j),x_dp(3,j)+bbR(3,2,j)],'b','LineWidth',1.2);
        platform(3) = plot3([x_dp(1,j)+bbR(1,3,j),x_dp(1,j)+bbR(1,4,j)],[x_dp(2,j)+bbR(2,3,j),x_dp(2,j)+bbR(2,4,j)],[x_dp(3,j)+bbR(3,3,j),x_dp(3,j)+bbR(3,4,j)],'b','LineWidth',1.2);
        platform(4) = plot3([x_dp(1,j)+bbR(1,1,j),x_dp(1,j)+bbR(1,4,j)],[x_dp(2,j)+bbR(2,1,j),x_dp(2,j)+bbR(2,4,j)],[x_dp(3,j)+bbR(3,1,j),x_dp(3,j)+bbR(3,4,j)],'b','LineWidth',1.2);
        % side shape
        platform(5) = plot3([x_dp(1,j)+bbR(1,1,j),x_dp(1,j)+bbR(1,5,j)],[x_dp(2,j)+bbR(2,1,j),x_dp(2,j)+bbR(2,5,j)],[x_dp(3,j)+bbR(3,1,j),x_dp(3,j)+bbR(3,5,j)],'b','LineWidth',1.2);
        platform(6) = plot3([x_dp(1,j)+bbR(1,6,j),x_dp(1,j)+bbR(1,2,j)],[x_dp(2,j)+bbR(2,6,j),x_dp(2,j)+bbR(2,2,j)],[x_dp(3,j)+bbR(3,6,j),x_dp(3,j)+bbR(3,2,j)],'b','LineWidth',1.2);
        platform(7) = plot3([x_dp(1,j)+bbR(1,3,j),x_dp(1,j)+bbR(1,7,j)],[x_dp(2,j)+bbR(2,3,j),x_dp(2,j)+bbR(2,7,j)],[x_dp(3,j)+bbR(3,3,j),x_dp(3,j)+bbR(3,7,j)],'b','LineWidth',1.2);
        platform(8) = plot3([x_dp(1,j)+bbR(1,8,j),x_dp(1,j)+bbR(1,4,j)],[x_dp(2,j)+bbR(2,8,j),x_dp(2,j)+bbR(2,4,j)],[x_dp(3,j)+bbR(3,8,j),x_dp(3,j)+bbR(3,4,j)],'b','LineWidth',1.2);
        % bottom shape
        platform(9) = plot3([x_dp(1,j)+bbR(1,6,j),x_dp(1,j)+bbR(1,5,j)],[x_dp(2,j)+bbR(2,6,j),x_dp(2,j)+bbR(2,5,j)],[x_dp(3,j)+bbR(3,6,j),x_dp(3,j)+bbR(3,5,j)],'b','LineWidth',1.2);
        platform(10)= plot3([x_dp(1,j)+bbR(1,6,j),x_dp(1,j)+bbR(1,7,j)],[x_dp(2,j)+bbR(2,6,j),x_dp(2,j)+bbR(2,7,j)],[x_dp(3,j)+bbR(3,6,j),x_dp(3,j)+bbR(3,7,j)],'b','LineWidth',1.2);
        platform(11)= plot3([x_dp(1,j)+bbR(1,8,j),x_dp(1,j)+bbR(1,7,j)],[x_dp(2,j)+bbR(2,8,j),x_dp(2,j)+bbR(2,7,j)],[x_dp(3,j)+bbR(3,8,j),x_dp(3,j)+bbR(3,7,j)],'b','LineWidth',1.2);
        platform(12)= plot3([x_dp(1,j)+bbR(1,8,j),x_dp(1,j)+bbR(1,5,j)],[x_dp(2,j)+bbR(2,8,j),x_dp(2,j)+bbR(2,5,j)],[x_dp(3,j)+bbR(3,8,j),x_dp(3,j)+bbR(3,5,j)],'b','LineWidth',1.2);
        %% mark part
        mark(1) = plot3([x_dp(1,j)+mmR(1,1,j),x_dp(1,j)+mmR(1,3,j)],[x_dp(2,j)+mmR(2,1,j),x_dp(2,j)+mmR(2,3,j)],[x_dp(3,j)+mmR(3,1,j),x_dp(3,j)+mmR(3,3,j)],'b','MarkerSize',8);
        mark(2) = plot3([x_dp(1,j)+mmR(1,2,j),x_dp(1,j)+mmR(1,4,j)],[x_dp(2,j)+mmR(2,2,j),x_dp(2,j)+mmR(2,4,j)],[x_dp(3,j)+mmR(3,2,j),x_dp(3,j)+mmR(3,4,j)],'b','MarkerSize',8); 
        mark(3) = plot3(x_dp(1,j)+mmR(1,1,j),x_dp(2,j)+mmR(2,1,j),x_dp(3,j)+mmR(3,1,j),'bo','MarkerSize',3); 
        mark(4) = plot3(x_dp(1,j)+mmR(1,2,j),x_dp(2,j)+mmR(2,2,j),x_dp(3,j)+mmR(3,2,j),'bo','MarkerSize',3); 
        mark(5) = plot3(x_dp(1,j)+mmR(1,3,j),x_dp(2,j)+mmR(2,3,j),x_dp(3,j)+mmR(3,3,j),'bo','MarkerSize',3); 
        mark(6) = plot3(x_dp(1,j)+mmR(1,4,j),x_dp(2,j)+mmR(2,4,j),x_dp(3,j)+mmR(3,4,j),'bo','MarkerSize',3); 
        mark(7) = plot3((x_dp(1,j)+mmR(1,1,j)+x_dp(1,j)+mmR(1,3,j))/2,(x_dp(2,j)+mmR(2,1,j)+x_dp(2,j)+mmR(2,3,j))/2,(x_dp(3,j)+mmR(3,1,j)+x_dp(3,j)+mmR(3,3,j))/2,'bo','MarkerSize',3); 
        %% Camera image part
        camera(1) = plot3([x_dp(1,j)+mmC(1,1,j),x_dp(1,j)+mmC(1,3,j)],[x_dp(2,j)+mmC(2,1,j),x_dp(2,j)+mmC(2,3,j)],[mmC(3,1,j)+x(3,j)-0.1,mmC(3,3,j)+x(3,j)-0.1],'r','MarkerSize',8);
        camera(2) = plot3([x_dp(1,j)+mmC(1,2,j),x_dp(1,j)+mmC(1,4,j)],[x_dp(2,j)+mmC(2,2,j),x_dp(2,j)+mmC(2,4,j)],[mmC(3,2,j)+x(3,j)-0.1,mmC(3,4,j)+x(3,j)-0.1],'r','MarkerSize',8); 
        camera(3) = plot3(x_dp(1,j)+mmC(1,1,j),x_dp(2,j)+mmC(2,1,j),mmC(3,1,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
        camera(4) = plot3(x_dp(1,j)+mmC(1,2,j),x_dp(2,j)+mmC(2,2,j),mmC(3,2,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
        camera(5) = plot3(x_dp(1,j)+mmC(1,3,j),x_dp(2,j)+mmC(2,3,j),mmC(3,3,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
        camera(6) = plot3(x_dp(1,j)+mmC(1,4,j),x_dp(2,j)+mmC(2,4,j),mmC(3,4,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
        camera(7) = plot3((x_dp(1,j)+mmC(1,1,j)+x_dp(1,j)+mmC(1,3,j))/2,(x_dp(2,j)+mmC(2,1,j)+x_dp(2,j)+mmC(2,3,j))/2,(mmC(3,1,j)+mmC(3,3,j)+x(3,j)-0.1+x(3,j)-0.1)/2,'ro','MarkerSize',3); 
        drawnow
        delete(p1);delete(p2);delete(p3);delete(platform);delete(mark);delete(camera)
    end
        %% quadrotor part
        p1=plot3(x(1,j),x(2,j),x(3,j),'g*');
        p2=plot3([x(1,j)+LLR(1,1,j),x(1,j)+LLR(1,3,j)],[x(2,j)+LLR(2,1,j),x(2,j)+LLR(2,3,j)],[x(3,j)+LLR(3,1,j),x(3,j)+LLR(3,3,j)],'k');
        p3=plot3([x(1,j)+LLR(1,2,j),x(1,j)+LLR(1,4,j)],[x(2,j)+LLR(2,2,j),x(2,j)+LLR(2,4,j)],[x(3,j)+LLR(3,2,j),x(3,j)+LLR(3,4,j)],'k');
        %% platform part
%         pp1=plot3(x_dp(1,j),x_dp(2,j),x_dp(3,j),'b*','MarkerSize',1);
        % upper shpae
        platform(1) = plot3([x_dp(1,j)+bbR(1,1,j),x_dp(1,j)+bbR(1,2,j)],[x_dp(2,j)+bbR(2,1,j),x_dp(2,j)+bbR(2,2,j)],[x_dp(3,j)+bbR(3,1,j),x_dp(3,j)+bbR(3,2,j)],'b','LineWidth',1.2);
        platform(2) = plot3([x_dp(1,j)+bbR(1,3,j),x_dp(1,j)+bbR(1,2,j)],[x_dp(2,j)+bbR(2,3,j),x_dp(2,j)+bbR(2,2,j)],[x_dp(3,j)+bbR(3,3,j),x_dp(3,j)+bbR(3,2,j)],'b','LineWidth',1.2);
        platform(3) = plot3([x_dp(1,j)+bbR(1,3,j),x_dp(1,j)+bbR(1,4,j)],[x_dp(2,j)+bbR(2,3,j),x_dp(2,j)+bbR(2,4,j)],[x_dp(3,j)+bbR(3,3,j),x_dp(3,j)+bbR(3,4,j)],'b','LineWidth',1.2);
        platform(4) = plot3([x_dp(1,j)+bbR(1,1,j),x_dp(1,j)+bbR(1,4,j)],[x_dp(2,j)+bbR(2,1,j),x_dp(2,j)+bbR(2,4,j)],[x_dp(3,j)+bbR(3,1,j),x_dp(3,j)+bbR(3,4,j)],'b','LineWidth',1.2);
        % side shape
        platform(5) = plot3([x_dp(1,j)+bbR(1,1,j),x_dp(1,j)+bbR(1,5,j)],[x_dp(2,j)+bbR(2,1,j),x_dp(2,j)+bbR(2,5,j)],[x_dp(3,j)+bbR(3,1,j),x_dp(3,j)+bbR(3,5,j)],'b','LineWidth',1.2);
        platform(6) = plot3([x_dp(1,j)+bbR(1,6,j),x_dp(1,j)+bbR(1,2,j)],[x_dp(2,j)+bbR(2,6,j),x_dp(2,j)+bbR(2,2,j)],[x_dp(3,j)+bbR(3,6,j),x_dp(3,j)+bbR(3,2,j)],'b','LineWidth',1.2);
        platform(7) = plot3([x_dp(1,j)+bbR(1,3,j),x_dp(1,j)+bbR(1,7,j)],[x_dp(2,j)+bbR(2,3,j),x_dp(2,j)+bbR(2,7,j)],[x_dp(3,j)+bbR(3,3,j),x_dp(3,j)+bbR(3,7,j)],'b','LineWidth',1.2);
        platform(8) = plot3([x_dp(1,j)+bbR(1,8,j),x_dp(1,j)+bbR(1,4,j)],[x_dp(2,j)+bbR(2,8,j),x_dp(2,j)+bbR(2,4,j)],[x_dp(3,j)+bbR(3,8,j),x_dp(3,j)+bbR(3,4,j)],'b','LineWidth',1.2);
        % bottom shape
        platform(9) = plot3([x_dp(1,j)+bbR(1,6,j),x_dp(1,j)+bbR(1,5,j)],[x_dp(2,j)+bbR(2,6,j),x_dp(2,j)+bbR(2,5,j)],[x_dp(3,j)+bbR(3,6,j),x_dp(3,j)+bbR(3,5,j)],'b','LineWidth',1.2);
        platform(10)= plot3([x_dp(1,j)+bbR(1,6,j),x_dp(1,j)+bbR(1,7,j)],[x_dp(2,j)+bbR(2,6,j),x_dp(2,j)+bbR(2,7,j)],[x_dp(3,j)+bbR(3,6,j),x_dp(3,j)+bbR(3,7,j)],'b','LineWidth',1.2);
        platform(11)= plot3([x_dp(1,j)+bbR(1,8,j),x_dp(1,j)+bbR(1,7,j)],[x_dp(2,j)+bbR(2,8,j),x_dp(2,j)+bbR(2,7,j)],[x_dp(3,j)+bbR(3,8,j),x_dp(3,j)+bbR(3,7,j)],'b','LineWidth',1.2);
        platform(12)= plot3([x_dp(1,j)+bbR(1,8,j),x_dp(1,j)+bbR(1,5,j)],[x_dp(2,j)+bbR(2,8,j),x_dp(2,j)+bbR(2,5,j)],[x_dp(3,j)+bbR(3,8,j),x_dp(3,j)+bbR(3,5,j)],'b','LineWidth',1.2);
        %% mark part
        mark(1) = plot3([x_dp(1,j)+mmR(1,1,j),x_dp(1,j)+mmR(1,3,j)],[x_dp(2,j)+mmR(2,1,j),x_dp(2,j)+mmR(2,3,j)],[x_dp(3,j)+mmR(3,1,j),x_dp(3,j)+mmR(3,3,j)],'b','MarkerSize',8);
        mark(2) = plot3([x_dp(1,j)+mmR(1,2,j),x_dp(1,j)+mmR(1,4,j)],[x_dp(2,j)+mmR(2,2,j),x_dp(2,j)+mmR(2,4,j)],[x_dp(3,j)+mmR(3,2,j),x_dp(3,j)+mmR(3,4,j)],'b','MarkerSize',8); 
        mark(3) = plot3(x_dp(1,j)+mmR(1,1,j),x_dp(2,j)+mmR(2,1,j),x_dp(3,j)+mmR(3,1,j),'bo','MarkerSize',3);  
        mark(4) = plot3(x_dp(1,j)+mmR(1,2,j),x_dp(2,j)+mmR(2,2,j),x_dp(3,j)+mmR(3,2,j),'bo','MarkerSize',3); 
        mark(5) = plot3(x_dp(1,j)+mmR(1,3,j),x_dp(2,j)+mmR(2,3,j),x_dp(3,j)+mmR(3,3,j),'bo','MarkerSize',3); 
        mark(6) = plot3(x_dp(1,j)+mmR(1,4,j),x_dp(2,j)+mmR(2,4,j),x_dp(3,j)+mmR(3,4,j),'bo','MarkerSize',3); 
        mark(7) = plot3((x_dp(1,j)+mmR(1,1,j)+x_dp(1,j)+mmR(1,3,j))/2,(x_dp(2,j)+mmR(2,1,j)+x_dp(2,j)+mmR(2,3,j))/2,(x_dp(3,j)+mmR(3,1,j)+x_dp(3,j)+mmR(3,3,j))/2,'bo','MarkerSize',3); 
        %% Camera image part
        camera(1) = plot3([x_dp(1,j)+mmC(1,1,j),x_dp(1,j)+mmC(1,3,j)],[x_dp(2,j)+mmC(2,1,j),x_dp(2,j)+mmC(2,3,j)],[mmC(3,1,j)+x(3,j)-0.1,mmC(3,3,j)+x(3,j)-0.1],'r','MarkerSize',8);
        camera(2) = plot3([x_dp(1,j)+mmC(1,2,j),x_dp(1,j)+mmC(1,4,j)],[x_dp(2,j)+mmC(2,2,j),x_dp(2,j)+mmC(2,4,j)],[mmC(3,2,j)+x(3,j)-0.1,mmC(3,4,j)+x(3,j)-0.1],'r','MarkerSize',8); 
        camera(3) = plot3(x_dp(1,j)+mmC(1,1,j),x_dp(2,j)+mmC(2,1,j),mmC(3,1,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
        camera(4) = plot3(x_dp(1,j)+mmC(1,2,j),x_dp(2,j)+mmC(2,2,j),mmC(3,2,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
        camera(5) = plot3(x_dp(1,j)+mmC(1,3,j),x_dp(2,j)+mmC(2,3,j),mmC(3,3,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
        camera(6) = plot3(x_dp(1,j)+mmC(1,4,j),x_dp(2,j)+mmC(2,4,j),mmC(3,4,j)+x(3,j)-0.1,'ro','MarkerSize',3); 
        camera(7) = plot3((x_dp(1,j)+mmC(1,1,j)+x_dp(1,j)+mmC(1,3,j))/2,(x_dp(2,j)+mmC(2,1,j)+x_dp(2,j)+mmC(2,3,j))/2,(mmC(3,1,j)+mmC(3,3,j)+x(3,j)-0.1+x(3,j)-0.1)/2,'ro','MarkerSize',3); 
end

%% RK4 method
function [pos,vel] = rk4_ode(X,m,R,U,g,Ig,dt)
    f1=ode(X,m,R,U,g,Ig);
    f2=ode(X+dt/2 *f1,m,R,U,g,Ig);
    f3=ode(X+dt/2 *f2,m,R,U,g,Ig);
    f4=ode(X+dt * f3,m,R,U,g,Ig);
    XX = (X + dt * (f1/6 + f2/3 +f3/3 + f4/6));
    pos = XX(1:6,1); vel = XX(7:12,1);
end

%% equations
function x_prime = ode(X,m,R,U,g,Ig)
   x_prime = [X(7:12); 
              ( -sin(X(4))*sin(X(6))-cos(X(4))*sin(X(5))*cos(X(6)) ) * U(1)/m;
              ( sin(X(4))*cos(X(6))-cos(X(4))*sin(X(5))*sin(X(6)) ) * U(1)/m;
              g - (cos(X(4))*cos(X(5))) * U(1)/m;
              U(2)/Ig(1,1) + (Ig(2,2)-Ig(3,3))*X(11)*X(12)/Ig(1,1);
              U(3)/Ig(2,2) + (Ig(3,3)-Ig(1,1))*X(10)*X(12)/Ig(2,2);
              U(4)/Ig(3,3) + (Ig(1,1)-Ig(2,2))*X(10)*X(11)/Ig(3,3)];
end
