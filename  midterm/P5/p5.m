%Tyler Matthews
%System Simluation Midterm P5
clc; close all; %Clear Console and Close Figures


%% PART A
N = 10000;
t = linspace(0,10,N);
u = ones(1,N);
stable_acc = 0.146; %Stable for predictor


fx1c = zeros(1,N);
fx2c = zeros(1,N);
fx3c = zeros(1,N);

x1p = zeros(1,N);
x2p = zeros(1,N);
x3p = zeros(1,N);

fx1p = zeros(1,N);
fx2p = zeros(1,N);
fx3p = zeros(1,N);

x1c = zeros(1,N);
x2c = zeros(1,N);
x3c = zeros(1,N);
y = zeros(1,N);

x1c(1) = 1;
x2c(1) = 1;
x3c(1) = 1;
x1p(1) = 1;
x2p(1) = 1;
x3p(1) = 1;

for k = 1:N-2
    %Predict
    x1p(k+2) = 1.45*x1c(k+1) - 0.45*x1c(k) + stable_acc * (1.27*fx1c(k+1) - 0.73*fx1c(k));
    x2p(k+2) = 1.45*x2c(k+1) - 0.45*x2c(k) + stable_acc * (1.27*fx2c(k+1) - 0.73*fx2c(k));
    x3p(k+2) = 1.45*x3c(k+1) - 0.45*x3c(k) + stable_acc * (1.27*fx3c(k+1) - 0.73*fx3c(k));
    
    fx1p(k+2) = -4.7*x1p(k+2)-1.55*x2p(k+2)-0.55*x3p(k+2)+u(k+2);
    fx2p(k+2) = 0.3*x1p(k+2)-2.75*x2p(k+2)-0.35*x3p(k+2);
    fx3p(k+2) = 1.1*x1p(k+2)+1.85*x2p(k+2)-2.55*x3p(k+2)-u(k+2);

    %Correct
    x1c(k+2) = 1.56*x1c(k+1) - 0.56*x1c(k) + stable_acc * (0.46*fx1p(k+2) + 0.29*fx1p(k+1) - 0.32*fx1p(k));
    x2c(k+2) = 1.56*x2c(k+1) - 0.56*x2c(k) + stable_acc * (0.46*fx2p(k+2) + 0.29*fx2p(k+1) - 0.32*fx2p(k));
    x3c(k+2) = 1.56*x3c(k+1) - 0.56*x3c(k) + stable_acc * (0.46*fx3p(k+2) + 0.29*fx3p(k+1) - 0.32*fx3p(k));
    
    fx1c(k+2) = -4.7*x1c(k+2)-1.55*x2c(k+2)-0.55*x3c(k+2)+u(k+2);
    fx2c(k+2) = 0.3*x1c(k+2)-2.75*x2c(k+2)-0.35*x3c(k+2);
    fx3c(k+2) = 1.1*x1c(k+2)+1.85*x2c(k+2)-2.55*x3c(k+2)-u(k+2);
    
 
    y(k) = 2*x1c(k+1)+x2c(k+1)+x3c(k+1);
end

figure
plot(t,y)
xlim([0 0.2])
title('Predictor Corrector Pair With T from P3')


%% PART B

N = 10000;
t = linspace(0,10,N);
u = ones(1,N);
stable_acc = 0.23; %stable for corrector


fx1c = zeros(1,N);
fx2c = zeros(1,N);
fx3c = zeros(1,N);

x1p = zeros(1,N);
x2p = zeros(1,N);
x3p = zeros(1,N);

fx1p = zeros(1,N);
fx2p = zeros(1,N);
fx3p = zeros(1,N);

x1c = zeros(1,N);
x2c = zeros(1,N);
x3c = zeros(1,N);
y = zeros(1,N);

x1c(1) = 1;
x2c(1) = 1;
x3c(1) = 1;
x1p(1) = 1;
x2p(1) = 1;
x3p(1) = 1;

for k = 1:N-2
    
    %Predict
    x1p(k+2) = 1.45*x1c(k+1) - 0.45*x1c(k) + stable_acc * (1.27*fx1c(k+1) - 0.73*fx1c(k));
    x2p(k+2) = 1.45*x2c(k+1) - 0.45*x2c(k) + stable_acc * (1.27*fx2c(k+1) - 0.73*fx2c(k));
    x3p(k+2) = 1.45*x3c(k+1) - 0.45*x3c(k) + stable_acc * (1.27*fx3c(k+1) - 0.73*fx3c(k));
    
    fx1p(k+2) = -4.7*x1p(k+2)-1.55*x2p(k+2)-0.55*x3p(k+2)+u(k+2);
    fx2p(k+2) = 0.3*x1p(k+2)-2.75*x2p(k+2)-0.35*x3p(k+2);
    fx3p(k+2) = 1.1*x1p(k+2)+1.85*x2p(k+2)-2.55*x3p(k+2)-u(k+2);

    %Correct
    x1c(k+2) = 1.56*x1c(k+1) - 0.56*x1c(k) + stable_acc * (0.46*fx1p(k+2) + 0.29*fx1p(k+1) - 0.32*fx1p(k));
    x2c(k+2) = 1.56*x2c(k+1) - 0.56*x2c(k) + stable_acc * (0.46*fx2p(k+2) + 0.29*fx2p(k+1) - 0.32*fx2p(k));
    x3c(k+2) = 1.56*x3c(k+1) - 0.56*x3c(k) + stable_acc * (0.46*fx3p(k+2) + 0.29*fx3p(k+1) - 0.32*fx3p(k));
    
    fx1c(k+2) = -4.7*x1c(k+2)-1.55*x2c(k+2)-0.55*x3c(k+2)+u(k+2);
    fx2c(k+2) = 0.3*x1c(k+2)-2.75*x2c(k+2)-0.35*x3c(k+2);
    fx3c(k+2) = 1.1*x1c(k+2)+1.85*x2c(k+2)-2.55*x3c(k+2)-u(k+2);
    
 
    y(k) = 2*x1c(k+1)+x2c(k+1)+x3c(k+1);
end

figure
plot(t,y)
xlim([0 0.2])
title('Predictor Corrector Pair With T from P4')

%% PART C
disp('The plots for Part A and Part B are similiar but obviously different')
disp('This is because the T value for PART A is both stable and accurate for the predictor and corrector')
disp('while the T value for PART B is stable and accurate for the correct, but stable and inaccurate for the predictor')
disp('This is because the stability region for the corrector is much larger and encompasses that of the corrector')
disp('Meaning that values that look to be stable and accurate for the corrector are not necessilarily acceptable for the predictor')

disp('This effect is exaggerated in Figure 3 where a stable and accurate value for the corretor is chosen that is unstable for the predictor')


N = 10000;
t = linspace(0,10,N);
u = ones(1,N);
stable_acc = 0.4280; %Stable for corrector, but not predictor


fx1c = zeros(1,N);
fx2c = zeros(1,N);
fx3c = zeros(1,N);

x1p = zeros(1,N);
x2p = zeros(1,N);
x3p = zeros(1,N);

fx1p = zeros(1,N);
fx2p = zeros(1,N);
fx3p = zeros(1,N);

x1c = zeros(1,N);
x2c = zeros(1,N);
x3c = zeros(1,N);
y = zeros(1,N);

x1c(1) = 1;
x2c(1) = 1;
x3c(1) = 1;
x1p(1) = 1;
x2p(1) = 1;
x3p(1) = 1;

for k = 1:N-2
    
    %Predict
    x1p(k+2) = 1.45*x1c(k+1) - 0.45*x1c(k) + stable_acc * (1.27*fx1c(k+1) - 0.73*fx1c(k));
    x2p(k+2) = 1.45*x2c(k+1) - 0.45*x2c(k) + stable_acc * (1.27*fx2c(k+1) - 0.73*fx2c(k));
    x3p(k+2) = 1.45*x3c(k+1) - 0.45*x3c(k) + stable_acc * (1.27*fx3c(k+1) - 0.73*fx3c(k));
    
    fx1p(k+2) = -4.7*x1p(k+2)-1.55*x2p(k+2)-0.55*x3p(k+2)+u(k+2);
    fx2p(k+2) = 0.3*x1p(k+2)-2.75*x2p(k+2)-0.35*x3p(k+2);
    fx3p(k+2) = 1.1*x1p(k+2)+1.85*x2p(k+2)-2.55*x3p(k+2)-u(k+2);

    %Correct
    x1c(k+2) = 1.56*x1c(k+1) - 0.56*x1c(k) + stable_acc * (0.46*fx1p(k+2) + 0.29*fx1p(k+1) - 0.32*fx1p(k));
    x2c(k+2) = 1.56*x2c(k+1) - 0.56*x2c(k) + stable_acc * (0.46*fx2p(k+2) + 0.29*fx2p(k+1) - 0.32*fx2p(k));
    x3c(k+2) = 1.56*x3c(k+1) - 0.56*x3c(k) + stable_acc * (0.46*fx3p(k+2) + 0.29*fx3p(k+1) - 0.32*fx3p(k));
    
    fx1c(k+2) = -4.7*x1c(k+2)-1.55*x2c(k+2)-0.55*x3c(k+2)+u(k+2);
    fx2c(k+2) = 0.3*x1c(k+2)-2.75*x2c(k+2)-0.35*x3c(k+2);
    fx3c(k+2) = 1.1*x1c(k+2)+1.85*x2c(k+2)-2.55*x3c(k+2)-u(k+2);
    
 
    y(k) = 2*x1c(k+1)+x2c(k+1)+x3c(k+1);
end

figure
plot(t,y)
xlim([0 0.2])
title('Stable Corrector & Unstable predictor')
