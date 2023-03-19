# Analysis-and-Control-of-SCARA-Manipulator
MATLAB programming simulate the, kinematics dynamics, kinematic and dynamic control of a SCARA Manipulator
%% Error
% The tau1, tau2, tau3 are not having 'g' term in them
% d3ddot term is not present in tau3
% one entire column of M matrix is 0
% inv(M) doesnt exist 


clear all;close all;clc
%% Robot Parameters
syms th1 th2 d3 real;
syms d1 d2 d4 a1 a2 real
% d1=0.15;
% d2=0.05;
% d4=0.2;
% a1=0.15;
% a2=0.15;

%% SCARA Robot Non Standard DH Table 3DOF
dh_table=[0 0 th1 d1;
    0 a1 th2 d2;
    pi a2 0 d3;
    0 0 0 d4];

%% Transformation Matrix
for i=1:4
    T(:,:,i)=[cos(dh_table(i,3)) -sin(dh_table(i,3)) 0 dh_table(i,2);
        sin(dh_table(i,3))*cos(dh_table(i,1)) cos(dh_table(i,3))*cos(dh_table(i,1)) -sin(dh_table(i,1)) -dh_table(i,4)*sin(dh_table(i,1));
        sin(dh_table(i,3))*sin(dh_table(i,1)) cos(dh_table(i,3))*sin(dh_table(i,1)) cos(dh_table(i,1)) dh_table(i,4)*cos(dh_table(i,1));
        0 0 0 1];
    
end

%% Rotation Matrices
R0_1 = T(1:3,1:3,1);
R1_2 = T(1:3,1:3,2);
R2_3 = T(1:3,1:3,3);
R3_4 = T(1:3,1:3,4);
%% Position Matrices
P0_1 = T(1:3,4,1);
P1_2 = T(1:3,4,2);
P2_3 = T(1:3,4,3);
P3_4 = T(1:3,4,4);
%% Forward Kinematics
T0_1 = T(:,:,1);
T0_2 = T(:,:,1) * T(:,:,2);
T0_3 = T(:,:,1) * T(:,:,2) * T(:,:,3);
T0_4 = T(:,:,1) * T(:,:,2) * T(:,:,3) * T(:,:,4);
R0_4 = T0_4(1:3,1:3);
P0_4 = T0_4(1:3,4);

x=T0_4(1,4);
y=T0_4(2,4);
z=T0_4(3,4);

x1=T0_1(1,4);
y1=T0_1(2,4);
z1=T0_1(3,4);

x2=T0_2(1,4);
y2=T0_2(2,4);
z2=T0_2(3,4);

x3=T0_3(1,4);
y3=T0_3(2,4);
z3=T0_3(3,4);


%% Velocity Kinematics
%% Angular velocity propagation
syms th1dot th2dot d3dot  real
w0 = [0;0;0]; %angulat velocity initial
w1 = R0_1'*(w0) + [0;0;th1dot];
w2 = R1_2'*(w1) + [0;0;th2dot];
w3 = R2_3'*(w2);
w4 = R3_4'*(w3);
%end effector angular velocities wrt base frame
w04 = R0_4* w4;

%% Linear velocity propogation
v0 =[0;0;0];
v1 = R0_1'*(v0 + cross(w0,P0_1));
v2 = R1_2'*(v1 + cross(w1,P1_2));
v3 = R2_3'*(v2 + cross(w2,P2_3)) + [0;0;d3dot];
v4 = R3_4'*(v3 + cross(w3,P3_4));

% end effector linear velocities wrt base frame
v04 = simplify(R0_1*R1_2*R2_3*R3_4 * v4);

%% Jacobian Matrix
% through velocity propogation 
J = equationsToMatrix(v04,[th1dot,th2dot,d3dot]);
% through partial dervatives
J1=[diff(x,th1) diff(x,th2) diff(x,d3);
   diff(y,th1) diff(y,th2) diff(y,d3);
   diff(z,th1) diff(z,th2) diff(z,d3)];
J1=simplify(J1);

%% Workspace analysis
detJ = simplify(det(J));
% for singularity

sol= solve([detJ],[th1,th2,d3]);
th1sol = vpa(sol.th1); th2sol = vpa(sol.th2); d3sol = vpa(sol.d3);



%% Dynamic model
%location of centre of mass of links
syms IC0_11 IC0_12 IC0_13 IC0_21 IC0_22 IC0_23 IC0_31 IC0_32 IC0_33 real
syms IC1_11 IC1_12 IC1_13 IC1_21 IC1_22 IC1_23 IC1_31 IC1_32 IC1_33 real
syms IC2_11 IC2_12 IC2_13 IC2_21 IC2_22 IC2_23 IC2_31 IC2_32 IC2_33 real
syms IC3_11 IC3_12 IC3_13 IC3_21 IC3_22 IC3_23 IC3_31 IC3_32 IC3_33 real

syms PC0x PC0y PC0z PC1x PC1y PC1z PC2x PC2y PC2z PC3x PC3y PC3z real

PC0_0 = [PC0x; PC0y; PC0z];
PC1_1 = [PC1x; PC1y; PC1z];
PC2_2 = [PC2x; PC2y; PC2z];
PC3_3 = [PC3x; PC3y; PC3z];

% PC0_0 = [0; 0.00; 135.774];
% PC1_1 = [224.441; 0.00; 25.415];
% PC2_2 = [210.376; 0.015; 78.37];
% PC3_3 = [0; 0.00; 265.00];

%Inertia matrix
IC0 = [IC0_11, IC0_12, IC0_13;
       IC0_21, IC0_22, IC0_23;
       IC0_31, IC0_32, IC0_33];
   
IC1 = [IC1_11, IC1_12, IC1_13;
       IC1_21, IC1_22, IC1_23;
       IC1_31, IC1_32, IC1_33];

IC2 = [IC2_11, IC2_12, IC2_13;
       IC2_21, IC2_22, IC2_23;
       IC2_31, IC2_32, IC2_33];

IC3 = [IC3_11, IC3_12, IC3_13;
       IC3_21, IC3_22, IC3_23;
       IC3_31, IC3_32, IC3_33];   

% IC0 = [3.640E+09, 0, -0.374;
%        0, 3.640E+09,0;
%        -0.374, 0, 2.488E+09];
%    
% IC1 = [5.521E+07, -214.772, -2.608E+07;
%        -214.772, 6.691E+08, 27.365;
%        -2.608E+07, 27.365, 7.040E+08];
% 
% IC2 = [1.446E+09, 2.438E+05, -3.325E+08;
%        2.438E+05, 5.189E+09, 35319.925;
%        -3.325E+08, 35319.925, 5.688E+09];
% 
% IC3 = [3.488E+08, 0, 0;
%        0, 3.488E+08, 0.00;
%        0, 0.00, 1.484E+07];


%% Inertial values and acceleration variables 
syms m1 m2 m3 m4 g th1ddot th2ddot d3ddot real

%% Angular acceleration vector 
alpha0 = [0;0;0];
alpha1 = R0_1'*(alpha0 + cross(w0,[0;0;th1dot])) + [0;0;th1ddot];
alpha2 = R1_2'*(alpha1 + cross(w1,[0;0;th2dot])) + [0;0;th2ddot];
alpha3 = R2_3'*(alpha2);
alpha4 = R3_4'*(alpha3);

%% Linear acceleration vector
a0 = [0;0;g];
a1 = R0_1'*(a0 +cross(alpha0,P0_1)+cross(w0, cross(w0,P0_1)));
a2 = R1_2'*(a1 +cross(alpha1,P1_2)+cross(w1, cross(w1,P1_2)));
a3 = R2_3'*(a2 +cross(alpha2,P2_3)+cross(w2, cross(w2,P2_3)))+[0;0;d3ddot] + 2*cross(w3,[0;0;d3ddot]);
a4 = R3_4'*(a3 +cross(alpha3,P3_4)+cross(w3, cross(w3,P3_4)));

%% Linear acceleration of centre of mass of links
ac1 = a1 + cross(alpha1, PC1_1) + cross(w1,cross(w1,PC1_1));
ac2 = a2 + cross(alpha2, PC2_2) + cross(w2,cross(w2,PC2_2));
ac3 = a3 + cross(alpha3, PC3_3) + cross(w3,cross(w3,PC3_3));

%% Inertial forces of the links
F1 = m1*ac1;
F2 = m2*ac2;
F3 = m3*ac3;


%% Inertial moments 
N1 = IC1*alpha0 + cross(w0,(IC1*w0));
N2 = IC2*alpha1 + cross(w1,(IC2*w1));
N3 = IC3*alpha2 + cross(w2,(IC3*w2));


%% End effector forces and moments 
f4 = [0;0;0];
n4 = [0;0;0];

%% Joint forces
f3 = R3_4*f4 + F3;
f2 = R2_3*f3 + F2;
f1 = R1_2*f2 + F1;
f0 = R0_1*f1;

%% Joint Moments
n3 = R3_4*n4 + N3 + cross(PC3_3,F3) + cross(P3_4,(R3_4*f4));
n2 = R2_3*n3 + N2 + cross(PC2_2,F2) + cross(P2_3,(R2_3*f3));
n1 = R1_2*n2 + N1 + cross(PC1_1,F1) + cross(P1_2,(R1_2*f2));
n0 = R0_1*n1 + cross(P0_1,(R0_1*f1));

%% Vector of Inputs
tau1 = simplify(n1(3));
tau2 = simplify(n2(3));
tau3 = simplify(f3(3));
