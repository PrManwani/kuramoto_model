%lop 5 
%finding Kc
clear;
tol = 0.00001;
N=80;
T=200;
tau=0.1;
r_final = zeros(1,1000);
flag_r = 0;
flag_t = 0;
r_cos = zeros(1,T);
r_sin = zeros(1,T);
theta = zeros(N,T);
theta_dot = zeros(N,T);
K = 1;
while(K<=1000 && (flag_t==0 || flag_r==0))    
    K
    w = random('Normal',0,1,1,N);
    a = random('Normal',0,1,1,N);
    
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
    error1 = std(theta_dot(T,(1:80)));
        if(error1<tol)
            Kc_t = K
            flag_t = flag_t + 1;
        end
    
    error2 = std(r(T-50:T))
    if(error2<tol)
        Kc_r = K
        flag_r = flag_r + 1;
    end
    r_final(K) = mean(r(T-50:T));
    K = K+1
end
    
