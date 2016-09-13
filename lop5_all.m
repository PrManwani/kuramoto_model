%lop 5 
%k Vs r graph
%finding critical value of K = Kc
%Using both r and theta_dot
clear;
flag_r = 0
flag_t = 0
tol = 0.000001;
N=80;
T=1000;
tau=0.1;
r_final = zeros(1,500);
r_cos = zeros(1,T);
r_sin = zeros(1,T);
theta = zeros(N,T);
theta_dot = zeros(N,T);
w = random('Normal',0,1,1,N);
a = random('Normal',0,1,1,N);
for K = 1:1:500
    
   
    for i=1:N
        theta(1,i)=a(i);
        theta_dot(1,i) = w(i);
    end
    for i =1:N
        for j = 1:N
            r_cos(1) = r_cos(1) + (1/N)*cos(theta(1,j)-theta(1,i));
            r_sin(1) = r_sin(1) + (1/N)*sin(theta(1,j)-theta(1,i));
        end
    end
    r_cos(1) = r_cos(1)/N;
    r_sin(1) = r_sin(1)/N;
    
    
    for t=2:T
        for i=1:N
        theta(t,i) = theta(t-1,i) + tau*theta_dot(t-1,i);
        theta(t,i) = mod(theta(t,i),2*pi);
            for j = 1:N
                theta_dot(t,i) = (w(i) + (K/N)*sin(theta(t,j)-theta(t,i)));         
            end
        end
        for i =1:N
            for j = 1:N
                r_cos(t) = r_cos(t) + (1/N)*cos(theta(t,j)-theta(t,i));
                r_sin(t) = r_sin(t) + (1/N)*sin(theta(t,j)-theta(t,i));
            end
        end
        r_cos(t) = r_cos(t)/N;
        r_sin(t) = r_sin(t)/N;
        
        r(t) = sqrt(r_cos(t)^2 + r_sin(t)^2);
        
        
    end
    theta_dot_k(K,(1:80)) = theta_dot(T,(1:80));
    theta_k(K,(1:80)) = theta(T,(1:80));    
    r_final(K) = mean(r(T-50:T));
    error1 = std(theta_dot(T,(1:80)));
    if(error1<tol && flag_t==0)
        Kc_t = K; 
        flag_t = flag_t + 1;
        
    end    
    error2 = std(r(T-50:T));
    if(error2<tol && flag_r ==0)
        Kc_r = K;
        flag_r = flag_r + 1;
    end
end
K = 1:1:500;


figure(1)
title('theta_dot vs K')
plot(K,theta_dot_k)
figure(2)
title('theta vs K')
plot(K,theta_k)
figure(3)
title('r vs K')
plot(K,r_final)
