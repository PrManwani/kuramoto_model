%plot theta_dot vs time for 8 different ks for different oscillators
clear;
N=80;
T=2400;
tau=0.1
w = random('Normal',pi/2,1,1,N);

K = [0 0.5 1 1.5 3 4 5 9]
k = zeros(8,10);
% for i=1:N
%         theta(1,i)=rand(;
% end
theta = rand(N)
rx=0;
ry=0;
phi(1)=0;
for i=1:N
    phi(1) = phi(1) + (1/N)*theta(1,i); 
    rx=rx+(1/N)*cos(theta(1,i)); 
    ry=ry+(1/N)*sin(theta(1,i)); 
end    

r(1) = sqrt(rx*rx + ry*ry);
r(2) = r(1);
phi(2)=phi(1);
for t=2:T
    rx=0;
    ry=0;
    phi(t+1)=0;
    for i=1:N
        k(i) = K(mod(i+7,8) + 1);
        theta_dot(t,i) = (w(i) + k(i)*r(t)*sin(phi(t)-theta(t-1,i)));
        theta(t,i) = theta(t-1,i) + tau*(w(i) + k(i)*r(t)*sin(phi(t)-theta(t-1,i)));
        theta(t,i) = mod(theta(t,i),2*pi);
        rx=rx+(1/N)*cos(theta(t,i));
        ry=ry+(1/N)*sin(theta(t,i)); 
        phi(t+1) = phi(t+1) + (1/N)*theta(t,i);
    end
    r(t+1) = sqrt(rx*rx + ry*ry);

    hold on
end
t = 1:1:T
j = 1
%no. of oscillators to show
m = 8
f = T/12%multiplication factor
figure(1)
title('theta_dot vs t')
subplot(3,2,1)
plot(t((1):(f*j)), theta_dot(t((1):(f*j)),(m:8:32)));
for j =2:6
subplot(3,2,j)
plot(t((f*(j-1)):(f*j)), theta_dot(t((f*(j-1)):(f*j)),(m:8:32)));
end
figure(2)
for j =7:12
subplot(3,2,j-6)
plot(t((f*(j-1)):(f*j)), theta_dot(t((f*(j-1)):(f*j)),(m:8:32)));
end