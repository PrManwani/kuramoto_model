%To note differenc between mean field coupling and actual model
%r is the mean field  coupling parameter
%K is the coupling constant
clear;
N=80;
T=1000;
tau=0.1;
phi(1)=0;
K = 290


w = random('Normal',0,1,1,N);
a = random('Normal',0,1,1,N);
r_cos = zeros(1,T);
r_sin = zeros(1,T);
for i=1:N
    theta(1,i)=a(i);
    theta_dot(1,i) = w(i);
    
end
for j = 1:N
    r_cos(1) = r_cos(1) + (1/N)*cos(theta(1,j)-theta(1,i));
    r_sin(1) = r_sin(1) + (1/N)*sin(theta(1,j)-theta(1,i));
end

for t=2:T
    for i=1:N
        theta(t,i) = theta(t-1,i) + tau*theta_dot(t-1,i);
        theta(t,i) = mod(theta(t,i),2*pi);
        theta_dot(t,i) = theta_dot(t-1,i);
        for j = 1:N
           theta_dot(t,i) = theta_dot(t,i) + (K/N)*sin(theta(t,j)-theta(t,i));         
        end
    end
    for j = 1:N
        r_cos(t) = r_cos(t) + (1/N)*cos(theta(t,j));
        r_sin(t) = r_sin(t) + (1/N)*sin(theta(t,j));
    end

    r(t) = sqrt(r_cos(t)^2 + r_sin(t)^2); 
    
  
end 
t = 1:1:T;
b = ((theta_dot(:,1)')); 
%figure(1)
%plot(t,b);
figure(2)
plot(t,r);