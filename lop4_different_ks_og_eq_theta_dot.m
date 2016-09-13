%3 different ks for all oscillators
%plot theta_dot vs time for 8 different ks for different oscillators

clear;
K = [0 50 100 150 200 250 300 350]
N=80;
T=1000;
tau=0.1

0;
w = random('Normal',0,1,1,N);
a = random('Normal',0,1,1,N);
r_cos = zeros(1,T);
r_sin = zeros(1,T);
for i=1:N
    theta(1,i)=a(i);
    theta_dot(1,i) = w(i);
    
end
for i =1:N
    for j = 1:N
        r_cos(1) = r_cos(1) + (1/N)*cos(theta(1,j)-theta(1,i));
        r_sin(1) = r_sin(1) + (1/N)*cos(theta(1,j)-theta(1,i));
    end
end    

for i =1:N
    for j = 1:N
        r_cos(1) = r_cos(1) + (1/N)*cos(theta(1,j)-theta(1,i));
        r_sin(1) = r_sin(1) + (1/N)*cos(theta(1,j)-theta(1,i));
    end
end    

for t=2:T
    for i=1:N
        k(i) = K(mod(i+7,8) + 1);
        theta(t,i) = theta(t-1,i) + tau*theta_dot(t-1,i);
        theta(t,i) = mod(theta(t,i),2*pi);
            for j = 1:N
                theta_dot(t,i) = (w(i) + (K/N)*sin(theta(t,j)-theta(t,i)));         
            end
    end
    for i =1:N
        for j = 1:N
        r_cos(t) = r_cos(t) + (1/N)*cos(theta(t,j)-theta(t,i));
        r_sin(t) = r_sin(t) + (1/N)*cos(theta(t,j)-theta(t,i));
        end
    end    
  
    r(t) = sqrt(r_cos(t)^2 + r_sin(t)^2); 

end
t = 1:1:T
j = 1
%no. of oscillators to show
m = 8
f = T/12%multiplication factor
l = 3:8:19
figure(1)
title('theta_dot vs t')
subplot(3,2,1)
plot(t((1):(f*j)), theta_dot(t((1):(f*j)),(1:8:17)));
for j =2:6
subplot(3,2,j)
plot(t((f*(j-1)):(f*j)), theta_dot(t((f*(j-1)):(f*j)),(1:8:17)));
end
figure(2)
for j =7:12
subplot(3,2,j-6)
plot(t((f*(j-1)):(f*j)), theta_dot(t((f*(j-1)):(f*j)),(1:8:17)));
end
