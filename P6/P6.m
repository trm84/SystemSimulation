%Tyler Matthews 2/23/2019
clc; %Clear console

%PART A
%Constants
a = 0.2;
b = 0.2;
c = 5.7;

for i=0 : 2
    %TIME STEP
    T = 0.0001;
    temp = 10 * 10^i;
    t = 0:T:temp;
    N = length(t);

    x = zeros(1,N);
    y = zeros(1,N);
    z = zeros(1,N);

    x(1) = 5;
    y(1) = 5;
    z(1) = 10;

    fx = (-1*y(1))-z(1);
    fy = x(1) + a*y(1);
    fz = b + z(1)*(x(1)-c);

    x(2) = x(1) + T*fx(1);
    y(2) = y(1) + T*fy(1);
    z(2) = z(1) + T*fz(1);

    for k=1 : N-1
        fx(k+1) = -1*(y(k+1)) - z(k+1);
        x(k+2)  = x(k+1) + (3/2)*T*fx(k+1) - (T/2)*fx(k);
        fy(k+1) = x(k+1) + a*y(k+1);
        y(k+2)  = y(k+1) + (3/2)*T*fy(k+1) - (T/2)*fy(k); 
        fz(k+1) = b + z(k+1)*(x(k+1)-c);
        z(k+2)  = z(k+1) + (3/2)*T*fz(k+1) - (T/2)*fz(k); 
    end
    figure;
    plot(x,y);
    xlabel('x');
    ylabel('y');
    TITLE = sprintf("Part A - t = %i", temp);
    title(TITLE);
end

%PART B
%Constants
a = 0.3;
b = 0.3;
c = 8.0;

%Timestep
T = 0.0001;
t = 0:T:250;
N = length(t);

x = zeros(1,N);
y = zeros(1,N);
z = zeros(1,N);

x(1) = 5;
y(1) = 5;
z(1) = 10;

fx = (-1*y(1))-z(1);
fy = x(1) + a*y(1);
fz = b + z(1)*(x(1)-c);

x(2) = x(1) + T*fx(1);
y(2) = y(1) + T*fy(1);
z(2) = z(1) + T*fz(1);

for k=1 : N-1
    fx(k+1) = -1*(y(k+1)) - z(k+1);
    x(k+2)  = x(k+1) + (3/2)*T*fx(k+1) - (T/2)*fx(k);
    fy(k+1) = x(k+1) + a*y(k+1);
    y(k+2)  = y(k+1) + (3/2)*T*fy(k+1) - (T/2)*fy(k); 
    fz(k+1) = b + z(k+1)*(x(k+1)-c);
    z(k+2)  = z(k+1) + (3/2)*T*fz(k+1) - (T/2)*fz(k); 
end
%2D Plot
figure;
plot(x,y);
xlabel('x');
ylabel('y');
TITLE = sprintf("Part B - 2D @ t = 250");
title(TITLE);

%3D Plot
figure;
plot3(x,y,z);
xlabel('x');
ylabel('y');
zlabel('z');
TITLE = sprintf("Part B - 3D @ t = 250");
title(TITLE);

%PART C
%Constants
a = 0.2;
b = 0.2;
c = 5.7;

%TIME STEP
for i=0 : 1
    T = 0.01;
    temp = 1000001 * 10^i;
    t = 0:T:4096;
    N = temp;

    x = zeros(1,N);
    y = zeros(1,N);
    z = zeros(1,N);

    x(1) = 5;
    y(1) = 5;
    z(1) = 10;

    fx = (-1*y(1))-z(1);
    fy = x(1) + a*y(1);
    fz = b + z(1)*(x(1)-c);

    x(2) = x(1) + T*fx(1);
    y(2) = y(1) + T*fy(1);
    z(2) = z(1) + T*fz(1);

    for k=1 : N-1
        fx(k+1) = -1*(y(k+1)) - z(k+1);
        x(k+2)  = x(k+1) + (3/2)*T*fx(k+1) - (T/2)*fx(k);
        fy(k+1) = x(k+1) + a*y(k+1);
        y(k+2)  = y(k+1) + (3/2)*T*fy(k+1) - (T/2)*fy(k); 
        fz(k+1) = b + z(k+1)*(x(k+1)-c);
        z(k+2)  = z(k+1) + (3/2)*T*fz(k+1) - (T/2)*fz(k); 
    end
    figure;
    plot(x(500000:N),y(500000:N));
    xlabel('x');
    ylabel('y');
    TITLE = sprintf("Part C - N = %i", temp);
    title(TITLE);
end
