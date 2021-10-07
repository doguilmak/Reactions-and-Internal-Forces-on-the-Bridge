clc
clear

E = 210e6 % kPa
A = 0.75 % m²
I = 0.141 % m4
q = 20 % kN/m
rad2deg = 180 / pi; % degree

ks1 = SpringElementStiffness(1e5) % kN/m
ks2 = SpringElementStiffness(2e5) % kN/m

%Lengths
L1 = 70 % m
L2 = PlaneFrameElementLength(0, 0, 60, 26) % m
L3 = PlaneFrameElementLength(0, 0, 40, 4) % m
L4 = PlaneFrameElementLength(0, 0, 40, -4) % m
L5 = PlaneFrameElementLength(0, 0, 60, -26) % m
L6 = 70 % m
L7 = 65 % m
L8 = 35 % m
L9 = 35 % m
L10 = 65 % m
L11 = PlaneFrameElementLength(0, 0, -5, 26) % m
L12 = 30 % m
L13 = PlaneFrameElementLength(0, 0, 5, 26) % m

%Degrees
theta1 = atan(26 / 60) * rad2deg;
theta2 = atan(4 / 40) * rad2deg;
theta3 = 360 - theta2;
theta4 = 360 - theta1;
theta5 = 270 + atan(5 / 26) * rad2deg;
theta6 = 270 - atan(5 / 26) * rad2deg;

% Building global stiffness matrix for the structure.
k1 = PlaneFrameElementStiffness(E, A, I, L1, 0);
k2 = PlaneFrameElementStiffness(E, A, I, L2, theta1);
k3 = PlaneFrameElementStiffness(E, A, I, L3, theta2);
k4 = PlaneFrameElementStiffness(E, A, I, L4, theta3);
k5 = PlaneFrameElementStiffness(E, A, I, L5, theta4);
k6 = PlaneFrameElementStiffness(E, A, I, L6, 0);
k7 = PlaneFrameElementStiffness(E, A, I, L7, 0);
k8 = PlaneFrameElementStiffness(E, A, I, L8, 0);
k9 = PlaneFrameElementStiffness(E, A, I, L9, 0);
k10 = PlaneFrameElementStiffness(E, A, I, L10, 180);
k11 = PlaneFrameElementStiffness(E, A, I, L11, theta5);
k12 = PlaneFrameElementStiffness(E, A, I, L12, 270);
k13 = PlaneFrameElementStiffness(E, A, I, L13, theta6);

%The springs have also DoF: nodes 11, 12, 13 and 14 - every node has one
%DoF, so total number of DoFs = 10 * 3 + 4 * 1 = 34
K = zeros(34, 34);

K = PlaneFrameAssemble(K, k1, 1, 2);
K = PlaneFrameAssemble(K, k2, 2, 3);
K = PlaneFrameAssemble(K, k3, 3, 4);
K = PlaneFrameAssemble(K, k4, 4, 5);
K = PlaneFrameAssemble(K, k5, 5, 6);
K = PlaneFrameAssemble(K, k6, 6, 7);
K = PlaneFrameAssemble(K, k7, 2, 8);
K = PlaneFrameAssemble(K, k8, 8, 9);
K = PlaneFrameAssemble(K, k9, 9, 10);
K = PlaneFrameAssemble(K, k10, 6, 10);
K = PlaneFrameAssemble(K, k11, 3, 8);
K = PlaneFrameAssemble(K, k12, 4, 9);
K = PlaneFrameAssemble(K, k13, 5, 10);

K = SpringAssemble(K, ks1, 2, 31);
K = SpringAssemble(K, ks2, 5, 32);
K = SpringAssemble(K, ks2, 17, 33);
K = SpringAssemble(K, ks1, 20, 34);

%The last arguments of the SpringAssemble() function are node numbers 
%only for systems that contain only spring elements (no other element types).
%Since spring nodes have only one degree of freedom, you need to provide 
%freedom degree numbers instead of node numbers to properly assemble the springs.


%Now you need to extract 30 columns and 30 rows from the global stiffness matrix.
%Hint: displacements in spring supports are known and equal to zero - 
%you don't need columns and rows that are the equations of these supports.

k = K(1:30,1:30);
f = [0; -700; -8167; 0; -1350; 1125; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -1350; -1125; 0; -700; 8167; 0; -1000; 5000; 0; -700; 0; 0; -1000; -5000];
u = k\f

% Starting to setting up global nodal displacement vector U.
% The calculated displacements for nodes 1-10 should be supplemented with four displacements in the spring supports.
% (i.e. nodes: 11 12 13 14, one displacement per node, which actually means one zero per node.)
U = [u;0;0;0;0];
% Calculading the global nodal force vector F.
F = K * U % Reactions in 1st, 2nd, 6th, 7th, 8th, 9th and 10th nodes.
%THE END
