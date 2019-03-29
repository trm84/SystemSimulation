%Tyler Matthews 2/13/19
%System Simulation Problem 05

clc; %Clear Console

%Constants
R1 = 500;
R2 = 1000;
R3 = 1000;
C1 = 4.7e-6;
C2 = 4.7e-6;
C3 = 4.7e-6;
L = 2;

%Matricies
A = [-1/(C1*R2), 1/(C1*R2), 0, 1/(C1); 1/(C2*R2), (-1/(C2*R2))+(-1/(C2*R3)), 1/(C2*R3), 0; 0, 1/(C3*R3), -1/(C3*R3), 0; -1/(L), 0, 0, (-1*R1)/(L)]
B = [0; 0; 0; 1/(L)]
C = [0, 0, 1, 0]
D = [0]

%Make transfer function from state space
[NUM, DEN] = ss2tf(A, B, C, D);
TF = tf(NUM, DEN)

%Poles of transfer function
POLES = pole(TF)

%Eigenvalues of A matrix
eigenvalues = eig(A)