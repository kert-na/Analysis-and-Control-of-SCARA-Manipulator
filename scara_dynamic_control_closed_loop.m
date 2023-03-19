%% Inverse dynamic simulation of a SCARA robot
clear all; close all; clc;

%% Simulation parameters
dt =0.1; %stepsize
ts = 15; %total simulationn time
t = 0:dt:ts; %time span
global a1 a2 d1 d2 d4 m1 m2 m3 m4 g

syms d3
%% Control parameter 
kp = 4; kd = 5;

%% System parameters
m1 =1; m2 =1; m3=1; m4=1;
d1=0.5;
d2=0.05;
d4=0.1;
a1=0.45;
a2=0.62; 
g=0*9.81 % gravity

%% Moment of innertia matrix
IC0_11 = 0; IC0_12 = 0; IC0_13 = 0;
IC0_21 = 0; IC0_22 = 0; IC0_23 = 0;
IC0_31 = 0; IC0_32 = 0; IC0_33 = 0;

IC1_11 = 0; IC1_12 = 0; IC1_13 = 0;
IC1_21 = 0; IC1_22 = 0; IC1_23 = 0;
IC1_31 = 0; IC1_32 = 0; IC1_33 = 0;
   
IC2_11 = 0; IC2_12 = 0; IC2_13 = 0;
IC2_21 = 0; IC2_22 = 0; IC2_23 = 0;
IC2_31 = 0; IC2_32 = 0; IC2_33 = 0;

IC3_11 = 0; IC3_12 = 0; IC3_31 = 0;
IC3_21 = 0; IC3_22 = 0; IC3_23 = 0;
IC3_31 = 0; IC3_32 = 0; IC3_33 = 0;



% IC0_11 = 3.640E+09; IC0_12 = 0; IC0_13 = -0.374;
% IC0_21 = 0; IC0_22 = 3.640E+09; IC0_23 = 0;
% IC0_31 = -0.374; IC0_32 = 0; IC0_33 = 2.488E+09;
% 
% IC1_11 = 5.521E+07; IC1_12 = -214.772; IC1_13 = -2.608E+07;
% IC1_21 = -214.772; IC1_22 = 6.691E+08; IC1_23 = 27.365;
% IC1_31 = -2.608E+07; IC1_32 = 27.365;  IC1_33 = 7.040E+08;
%    
% IC2_11 = 1.446E+09; IC2_12 = 2.438E+05; IC2_13 = -3.325E+08;
% IC2_21 = 2.438E+05; IC2_22 = 5.189E+09; IC2_23 = 35319.925;
% IC2_31 = -3.325E+08; IC2_32 = 35319.925; IC2_33 = 5.688E+09;
% 
% IC3_11 = 3.488E+08; IC3_12 = 0; IC3_31 = 0;
% IC3_21 = 0; IC3_22 = 3.488E+08; IC3_23 = 0;
% IC3_31 = 0; IC3_32 = 0; IC3_33 = 1.484E+07;


%% COM values
PC0x = 0; PC0y = 0; PC0z = d1;
PC1x = a1; PC1y = 0; PC1z = d2;
PC2x = a2; PC2y = 0; PC1z = -d3;
PC3x = 0; PC3y = 0; PC3z = d4;

% PC0x = 0; PC0y = 0; PC0z = 135.774;
% PC1x = 0.224441; PC1y = 0; PC1z = 25.415;
% PC2x = 210.376; PC2y = 0.015; PC1z = 78.37;
% PC3x = 0; PC3y = 0; PC3z = 265.00;


b1=0.5; b2=00.5; b3=1;
c1 =1; c2 = 1; c3 = 2;
g = -9.81 %gravity

%% Initial conditions
q_c=[0;0;0];
qdot_c = [0;0;0];

%% Numerical integration starts here
for i=1:length(t)
    %%Desired value
%     mu_d (:,i) = [2+1*sin(0.2*t(i));
%                   2-1*cos(0.2*t(i));
%                   d1 + d2 - 0.2*t(i) - d4];
    mu_d=[0.3;0.2;0.4];
%     mu_dot_d(:,i) = [0.2*cos(0.2*t(i));
%                      0.2*sin(0.2*t(i));
%                      -0.2];
    mu_dot_d =[0;0;0];
%     mu_ddot_d(:,i) = [-0.2*0.2*sin(0.2*t(i));
%                       0.2*0.2*cos(0.2*t(i));
%                       0];
    mu_ddot_d=[0;0;0];
                  
%% Desired joint positions and velociteis
    q_d(2,i) = acos((mu_d(1)^2 + mu_d(2)^2 - a1^2 - a2^2)/(2*a1*a2));
    q_d(1,i) = (a2*(mu_d(2)*cos(q_d(2,i))-mu_d(1)*sin(q_d(2,i)))+mu_d(2)*a1)/(a2*(mu_d(2)*sin(q_d(2,i))+mu_d(1)*cos(q_d(2,i)))+mu_d(1)*a1);
    q_d(3,i) = d1+d2-d4-mu_d(3);
    
    Jd = [- a2*sin(q_d(1,i) + q_d(2,i)) - a1*sin(q_d(1,i)), -a2*sin(q_d(1,i) + q_d(2,i)),  0;
              a2*cos(q_d(1,i) + q_d(2,i)) + a1*cos(q_d(1,i)),  a2*cos(q_d(1,i) + q_d(2,i)),  0;
              0,                  0, -1];
    qdot_d(:,i)=(Jd)\mu_dot_d;
%% Initial conditions
    q(:,1) = q_d(:,1);
    qdot(:,1)=qdot_d(:,1);
    
    mu(1,i)=a2*cos(q_c(1,i) + q_c(2,i)) + a1*cos(q_c(1,i));
    mu(2,i)=a2*sin(q_c(1,i) + q_c(2,i)) + a1*sin(q_c(1,i));
    mu(3,i)=d1 + d2 - q_c(3,i) - d4;
    
    J = [- a2*sin(q_c(1,i) + q_c(2,i)) - a1*sin(q_c(1,i)), -a2*sin(q_c(1,i) + q_c(2,i)),  0;
           a2*cos(q_c(1,i) + q_c(2,i)) + a1*cos(q_c(1,i)),  a2*cos(q_c(1,i) + q_c(2,i)),  0;
           0,                  0, -1];
    
%% Desired Dynamic terms
    Md= [IC2_33 - IC3_33 + PC1x^2*m1 + PC1y^2*m1 + PC2x^2*m2 + PC2y^2*m2 + PC3x^2*m3 + PC3y^2*m3 + a1^2*m2 + a1^2*m3 + a2^2*m3 + 2*PC3x*a2*m3 + 2*PC2x*a1*m2*cos(q_d(2,i)) + 2*PC3x*a1*m3*cos(q_d(2,i)) - 2*PC2y*a1*m2*sin(q_d(2,i)) + 2*PC3y*a1*m3*sin(q_d(2,i)) + 2*a1*a2*m3*cos(q_d(2,i)), m2*PC2x^2 + a1*m2*cos(q_d(2,i))*PC2x + m2*PC2y^2 - a1*m2*sin(q_d(2,i))*PC2y + m3*PC3x^2 + 2*m3*PC3x*a2 + a1*m3*cos(q_d(2,i))*PC3x + m3*PC3y^2 + a1*m3*sin(q_d(2,i))*PC3y + m3*a2^2 + a1*m3*cos(q_d(2,i))*a2 - IC3_33,  0;
                                                      m2*PC2x^2 + a1*m2*cos(q_d(2,i))*PC2x + m2*PC2y^2 - a1*m2*sin(q_d(2,i))*PC2y + m3*PC3x^2 + 2*m3*PC3x*a2 + a1*m3*cos(q_d(2,i))*PC3x + m3*PC3y^2 + a1*m3*sin(q_d(2,i))*PC3y + m3*a2^2 + a1*m3*cos(q_d(2,i))*a2 + IC2_33 - IC3_33,                                                                                                             m2*PC2x^2 + m2*PC2y^2 + m3*PC3x^2 + 2*m3*PC3x*a2 + m3*PC3y^2 + m3*a2^2 - IC3_33,  0;
                                                                                                                                                                                                                                                         0,                                                                                                                                                                                           0, m3;];
    OE_d = [-a1*qdot_d(2,i)*(2*qdot_d(1,i) + qdot_d(2,i))*(PC2x*m2*sin(q_d(2,i)) + PC3x*m3*sin(q_d(2,i)) + a2*m3*sin(q_d(2,i)) + PC2y*m2*cos(q_d(2,i)) - PC3y*m3*cos(q_d(2,i)));
                   a1*qdot_d(1,i)^2*(PC2x*m2*sin(q_d(2,i)) + PC3x*m3*sin(q_d(2,i)) + a2*m3*sin(q_d(2,i)) + PC2y*m2*cos(q_d(2,i)) - PC3y*m3*cos(q_d(2,i)));
                                                                                                                          0];
    GE_d = [0;0;-g*m3];
    
    FE = [b1*qdot_d(1,i)+c1*sign(qdot_d(1,i)); b2*qdot_d(2,i)+c2*sign(qdot_d(2,i)); b3*qdot_d(3,i)+c3*sign(qdot_d(3,i));];
    
    Jdot_d  = [- a2*cos(q_d(1,i) + q_d(2,i)) - a1*cos(q_d(1,i)), -a2*cos(q_d(1,i) + q_d(2,i)), 0;
               - a2*sin(q_d(1,i) + q_d(2,i)) - a1*sin(q_d(1,i)), -a2*sin(q_d(1,i) + q_d(2,i)), 0;
                                               0,                  0, 0;];
                                                         
 
%% ACtual dynamic terms
    M = [IC2_33 - IC3_33 + PC1x^2*m1 + PC1y^2*m1 + PC2x^2*m2 + PC2y^2*m2 + PC3x^2*m3 + PC3y^2*m3 + a1^2*m2 + a1^2*m3 + a2^2*m3 + 2*PC3x*a2*m3 + 2*PC2x*a1*m2*cos(q_c(2,i)) + 2*PC3x*a1*m3*cos(q_c(2,i)) - 2*PC2y*a1*m2*sin(q_c(2,i)) + 2*PC3y*a1*m3*sin(q_c(2,i)) + 2*a1*a2*m3*cos(q_c(2,i)), m2*PC2x^2 + a1*m2*cos(q_c(2,i))*PC2x + m2*PC2y^2 - a1*m2*sin(q_c(2,i))*PC2y + m3*PC3x^2 + 2*m3*PC3x*a2 + a1*m3*cos(q_c(2,i))*PC3x + m3*PC3y^2 + a1*m3*sin(q_c(2,i))*PC3y + m3*a2^2 + a1*m3*cos(q_c(2,i))*a2 - IC3_33,  0;
                                                      m2*PC2x^2 + a1*m2*cos(q_c(2,i))*PC2x + m2*PC2y^2 - a1*m2*sin(q_c(2,i))*PC2y + m3*PC3x^2 + 2*m3*PC3x*a2 + a1*m3*cos(q_c(2,i))*PC3x + m3*PC3y^2 + a1*m3*sin(q_c(2,i))*PC3y + m3*a2^2 + a1*m3*cos(q_c(2,i))*a2 + IC2_33 - IC3_33,                                                                                                             m2*PC2x^2 + m2*PC2y^2 + m3*PC3x^2 + 2*m3*PC3x*a2 + m3*PC3y^2 + m3*a2^2 - IC3_33,  0;
                                                                                                                                                                                                                                                         0,                                                                                                                                                                                           0, m3;];
    
    OE = [-a1*qdot_c(2,i)*(2*qdot_c(1,i) + qdot_c(2,i))*(PC2x*m2*sin(q_c(2,i)) + PC3x*m3*sin(q_c(2,i)) + a2*m3*sin(q_c(2,i)) + PC2y*m2*cos(q_c(2,i)) - PC3y*m3*cos(q_c(2,i)));
                   a1*qdot_c(1,i)^2*(PC2x*m2*sin(q_c(2,i)) + PC3x*m3*sin(q_c(2,i)) + a2*m3*sin(q_c(2,i)) + PC2y*m2*cos(q_c(2,i)) - PC3y*m3*cos(q_c(2,i)));
                                                                                                                          0];
                                                                                                                                                                                                                                                     
    GE = [0;0;-g*m3];
    
    FE = [b1*qdot_c(1,i)+c1*sign(qdot_c(1,i)); b2*qdot_c(2,i)+c2*sign(qdot_c(2,i)); b3*qdot_c(3,i)+c3*sign(qdot_c(3,i));];
 
    
    Jdot  = [- a2*cos(q_c(1,i) + q_c(2,i)) - a1*cos(q_c(1,i)), -a2*cos(q_c(1,i) + q_c(2,i)), 0;
             - a2*sin(q_c(1,i) + q_c(2,i)) - a1*sin(q_c(1,i)), -a2*sin(q_c(1,i) + q_c(2,i)), 0;
                                               0,                  0, 0;];

                                           
    %Inertia matrix in task space
    M_mu_d = inv(Jd')*Md*inv(Jd);

    %other effects in task space
    n_mu_d = inv(Jd')*(OE_d+GE_d-Md*inv(Jd)*Jdot_d*qdot_d(:,i));
    
   
    %% Errors
    mu_tilda(:,i)=mu_d-mu(:,i);
    mu_dot_tilda(:,i)=mu_dot_d-J*qdot_c(:,i);
    
    
    %% Computed Torque Control
    %input vector in task-space
    F_c(:,i)=M_mu_d*(mu_ddot_d+kp*mu_tilda(:,i)+kd*mu_dot_tilda(:,i)) + n_mu_d;
    
    %inpput vector in joint space
    tau_c(:,i) = Jd'*F_c(:,i);
    
    %% Inverse Dynamics
    %acc vector
    qddot_c(:,i)=inv(M)*(tau_c(:,i)-GE-OE);
    
    %vel vector
    qdot_c(:,i+1)=qdot_c(:,i)+qddot_c(:,i)*dt;
    
    %pos vect
    q_c(:,i+1) = q_c(:,i)+qdot_c(:,i)*dt + 1/2*qddot_c(:,i)*dt^2;
                                           
end

%% Animation figgu
for i=1:10:length(t)
    x0 = 0; y0 = 0; z0 = 0;

    x1 = x0; y1 = y0; z1 = z0 + d1;

    x2 = x1 + a1*cos(q_c(1,i)); y2 = y1 + a1*sin(q_c(1,i)); z2 = z1;

    x3 = x2; y3= y2; z3 = z2 + d2;

    x4 = x3 + a2*cos(q_c(1,i)+q_c(2,i)); y4 = y3 + a2*sin(q_c(1,i)+q_c(2,i)); z4 = z3;

    x5 = x4; y5 = y4; z5 = z4 - q_c(3,i);

    x6 = x5; y6 = y5; z6 = z5 - d4;


    plot3([x0,x1,x2,x3,x4,x5,x6],[y0,y1,y2,y3,y4,y5,y6],[z0,z1,z2,z3,z4,z5,z6],'linewidth',3)
    grid on, set(gca,'fontsize',20)
    axis([-(a1+a2)-0.1,(a1+a2)+0.1,-(a1+a2)-0.1,(a1+a2)+0.1,-0.1,(d1+d2)+0.1])
    hold off,pause(1)
end