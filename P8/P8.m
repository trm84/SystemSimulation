%Tyler Matthews
%System Simulation P8

clc; close all; %Clear console and close figures
Nt=21;
Nr=12;

theta=linspace(0,2*pi,1001); 
rho=linspace(0.984,1,1001); 
tvec=linspace(0,2*pi,Nt);
rvec=linspace(0.5,1,Nr);

figure;

hold on
 plot(rho*0,4*rho-2,'k')
 plot(4*rho-3,rho*0,'k')
hold off

for k=1:length(rvec)
 z=0.984*exp(1i*theta);
 w=(12*z.^3-12*z.^2)./(23*z.^2-16*z+5);
 
 hold on
 plot(real(w), imag(w))
 hold off
end

for k=1:length(tvec)-1
 z=rho*exp(1i*tvec(k));
 w=(12*z.^3-12*z.^2)./(23*z.^2-16*z+5);
 
 hold on
 plot(real(w), imag(w))
 hold off
end

title('AB3 Stability Region')
axis([-0.7 0.7 -1 1])
xlabel('Real')
ylabel('Imaginary')
grid on


%t95 = linspace(-1.2052,1.2052,1001);
z95 = 0.95*exp(1i*theta);
w95 = (12*z95.^3-12*z95.^2)./(23*z95.^2-16*z95+5);

%t925 = linspace(-1.4288,1.4288,1001);
z925 = 0.925*exp(1i*theta);
w925 = (12*z925.^3-12*z925.^2)./(23*z925.^2-16*z925+5);

%t90 = linspace(-1.5844,1.5844,1001);
z90 = 0.90*exp(1i*theta);
w90 = (12*z90.^3-12*z90.^2)./(23*z90.^2-16*z90+5);

%t875 = linspace(-1.7249,1.7249,1001);
z875 = 0.875*exp(1i*theta);
w875 = (12*z875.^3-12*z875.^2)./(23*z875.^2-16*z875+5);

%t85 = linspace(-1.8617,1.8617,1001);
z85 = 0.85*exp(1i*theta);
w85 = (12*z85.^3-12*z85.^2)./(23*z85.^2-16*z85+5);

%t825 = linspace(-2.0006,2.0006,1001);
z825 = 0.825*exp(1i*theta);
w825 = (12*z825.^3-12*z825.^2)./(23*z825.^2-16*z825+5);

%t80 = linspace(-2.1452,2.1452,1001);
z80 = 0.80*exp(1i*theta);
w80 = (12*z80.^3-12*z80.^2)./(23*z80.^2-16*z80+5);

%t775 = linspace(-2.2992,2.2992,1001);
z775 = 0.775*exp(1i*theta);
w775 = (12*z775.^3-12*z775.^2)./(23*z775.^2-16*z775+5);

%t75 = linspace(-2.4670,2.4670,1001);
z75 = 0.75*exp(1i*theta);
w75 = (12*z75.^3-12*z75.^2)./(23*z75.^2-16*z75+5);

%t725 = linspace(-2.6578,2.6578,1001);
z725 = 0.725*exp(1i*theta);
w725 = (12*z725.^3-12*z725.^2)./(23*z725.^2-16*z725+5);

Nw = 101;
zrp=zeros(Nt,Nw);
wrp=zeros(Nt,Nw);

figure;
hold on
 plot(real(w95), imag(w95), 'r')
 plot(real(w925), imag(w925), 'r')
 plot(real(w90), imag(w90), 'r')
 plot(real(w875), imag(w875), 'r')
 plot(real(w85), imag(w85), 'r')
 plot(real(w825), imag(w825), 'r')
 plot(real(w80), imag(w80), 'r')
 plot(real(w775), imag(w775), 'r')
 plot(real(w75), imag(w75), 'r')
 plot(real(w725), imag(w725), 'r')
 
 for m=1:Nt
 plot(real(wrp(m,:)),imag(wrp(m,:)))
 end
hold off

title('AB3 Stability : 0.725 to 0.95')
xlabel('Real')
ylabel('Imaginary')
grid on