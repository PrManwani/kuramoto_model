%3 different ks for all oscillators
%plotted theta vs t for all oscillators
clear;
N=80;
T=1000;
Klow = 0;
Kup = 5;
Kstep = 0.1;
totsteps = round((Kup-Klow)/Kstep);
Ktot(1) = Klow;

tau=0.1
 
%w = random('Normal',0,1,1,N);
%w = trnd(1,1,N);
w = rand(1,N);


  for m = 1:3
    K(1) = 0.01;
    K(2) = 2.2;
    K(3) = 4;
    
    
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
            
            theta_dot(t,i) = (w(i) + K(m)*r(t)*sin(phi(t)-theta(t-1,i)));
            theta(t,i) = theta(t-1,i) + tau*(w(i) + K(m)*r(t)*sin(phi(t)-theta(t-1,i)));
            theta(t,i) = mod(theta(t,i),2*pi);
            rx=rx+(1/N)*cos(theta(t,i));
            ry=ry+(1/N)*sin(theta(t,i)); 
            phi(t+1) = phi(t+1) + (1/N)*theta(t,i);
        end
        
        r(t+1) = sqrt(rx*rx + ry*ry); 
        figure(m)
        plot(t, theta(t,(1:5)), '+');
        hold on
        
    end
    
  end
    
    
    
   