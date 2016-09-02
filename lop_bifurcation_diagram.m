%Bifurcation diagram for all oscillators having same k in one run.
%Run over all k and plotted coupled frequency with k and angle with k
%plotted using frequency as a normal distribution
clear;
N=80;
T=200;
Klow = 0;
Kup = 5;
Kstep = 0.1;
totsteps = round((Kup-Klow)/Kstep);
Ktot(1) = Klow;
for i=2:totsteps;
    Ktot(i) =  Ktot(i-1)+Kstep;
end
tau=0.1
 
w = random('Normal',0,1,1,N);
%w = trnd(0.1,1,N);
%w = rand(1,N)



for Steps=1:totsteps
    
    K=Klow+ Steps*Kstep;
    
    
    for i=1:N
        theta(1,i)=w(i);
    end

    
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
            
            theta_dot(t,i) = (w(i) + K*r(t)*sin(phi(t)-theta(t-1,i)));
            theta(t,i) = theta(t-1,i) + tau*(w(i) + K*r(t)*sin(phi(t)-theta(t-1,i)));
            theta(t,i) = mod(theta(t,i),2*pi);
            rx=rx+(1/N)*cos(theta(t,i));
            ry=ry+(1/N)*sin(theta(t,i)); 
            phi(t+1) = phi(t+1) + (1/N)*theta(t,i);
        end
        
        r(t+1) = sqrt(rx*rx + ry*ry); 

    end
    figure(2)
    plot(K,theta_dot(200,(1:80)), '+');
    axis([Klow Kup -30 30])
    hold on
    
    figure(1)
    plot(K,theta(200,(1:80)), '+');
    axis([Klow Kup -8 8])
    hold on
    
    
    end
