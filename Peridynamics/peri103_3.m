function [R] = peri103_3(N,B,I,delta_x,M,rho,L0,a,x,v,E1,E2,E3,c1,c2,c3,gamma_M,gamma_m)

V=a(1:(N));

residual = 0;
%%
C3 = sqrt(E3/rho);
s_dot = M*C3;

%%
for j = (B+1):(B)+I % Interior nodes 
    partBsum = 0;
    partCsum = 0;
    partA = 0;
    for i= 1:N % Goes over all the nodes
        partB = 0;
        partC = 0;
        if i-j==0 % don't write as x(i)-x(j) cause of numerical errors
            %fprintf('0');
            partB = 0;
            partC = 0;
        else

            dudx = (V(i)-V(j))/(x(i)-x(j));% Look into this
            if dudx <= gamma_M
                F = (E1*dudx + c1);
            elseif dudx > gamma_M && dudx <= gamma_m
                F = (E2*dudx + c2);
            else
                F = (E3*dudx + c3);
            end
            partB = (L0^2/E3)*(4/sqrt(pi()))*((x(i) - x(j))/(L0^3))*(exp(-((x(i)- x(j))/L0)^2))*F*(delta_x/L0);
            
            %% adding in the viscous term
            if i==1 % V(i-1) = V(i)
                partC_partial = (1/(x(i) - x(j)))*(((V(i+1)-V(i))/(delta_x)) ...
                    - ((V(j+1)-V(j-1))/(2*delta_x)))*((x(i) - x(j))/(L0^3))*(exp(-((x(i)-x(j))/L0)^2));
            elseif i==N % V(i+1) = V(i)
                partC_partial = (1/(x(i) - x(j)))*(((V(i)-V(i-1))/(delta_x)) ...
                    - ((V(j+1)-V(j-1))/(2*delta_x)))*((x(i) - x(j))/(L0^3))*(exp(-((x(i)-x(j))/L0)^2));
            else
                partC_partial = (1/(x(i) - x(j)))*(((V(i+1)-V(i-1))/(2*delta_x)) ...
                    - ((V(j+1)-V(j-1))/(2*delta_x)))*((x(i) - x(j))/(L0^3))*(exp(-((x(i)-x(j))/L0)^2));
            end
            partC = (-v*s_dot/C3)*(L0^3/E3)*partC_partial*(delta_x/L0);
        end
        
        %% Computing sum from i=1 to i=N
        partBsum = partBsum + partB;
        partCsum = partCsum + partC;
    end   
    
    partA = (L0/C3^2)*(s_dot^2)*(1/delta_x^2)*(V(j+1)-2*V(j)+V(j-1));
    
    residual = residual +( partA - partBsum - partCsum )^2;
end
R = (residual);

