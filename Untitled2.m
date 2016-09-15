%lop 5 
%finding Kc
clear;
tol = 0.0001;
temp1 = 1e+04
temp2 = 1e+04
N=80;
T=4000;
tau=0.001;

flag_r = 0;
flag_t = 0;


K = 10;

while(K<=1000)
    actual_k = K/10;
    theta = zeros(T,N);
    theta_dot = zeros(T,N);
    r_cos = zeros(1,T);
    r_sin = zeros(1,T);
    r = zeros(1,T);
    w = random('Normal',0,1,1,N);
    a = random('Normal',0,1,1,N);
    
    for i=1:N
        theta(1,i)=a(i);
        theta_dot(1,i) = w(i);
    end
    for j = 1:N
        r_cos(1) = r_cos(1) + (1/N)*cos(theta(1,j));
        r_sin(1) = r_sin(1) + (1/N)*sin(theta(1,j));
    end
    r(1) = sqrt(r_cos(1)^2 + r_sin(1)^2); 
    for t=2:T
        for i=1:N
        theta(t,i) = theta(t-1,i) + tau*theta_dot(t-1,i);
        theta(t,i) = mod(theta(t,i),2*pi);
        theta_dot(t,i) = w(i);
            for j = 1:N
                theta_dot(t,i) = theta_dot(t,i) + (actual_k/N)*sin(theta(t,j)-theta(t,i));         
            end
        end
        for j = 1:N
            r_cos(t) = r_cos(t) + (1/N)*cos(theta(t,j));
            r_sin(t) = r_sin(t) + (1/N)*sin(theta(t,j));
        end                
        r(t) = sqrt(r_cos(t)^2 + r_sin(t)^2);
        
    end
    error1 = std(theta_dot(T,(1:N)));
    if(error1<temp1)
        error11 = error1;
        temp1 = error1;
    end    
    if(error1<tol)
            Kc_t = actual_k;
            flag_t = flag_t + 1;
    end
    
    error2 = std(r(T-50:T));
    if(error2<temp2)
        error22 = error2;
        temp2 = error2;
    end
    if(error2<tol)
        Kc_r = actual_k;
        flag_r = flag_r + 1;
    end
    r_final((K/10)) = mean(r(T-50:T));
    K = K+10
end
K = 10:10:1000
plot(K,r_final)