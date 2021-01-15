
%Lily's code to get eigen values

L = 12; %in
rho = 0.0002505 %lb-sec^2/in^4
E = 10175000 %lb/in^2

Mt = 1.131*rho
St = 0.5655*rho
It = 23.124*rho

A = 1*(1/8); %wh

Izz = (1/8)^3/12 %wh^3/12

cM2 = rho*A*L/100800
cK2 = 4*E*Izz/L^3

%master mass equation
M2 = cM2*[19272 1458*L 5298 -642*L 0 0; 1458*L 172*L^2 642*L -73*L^2 0 0;...
    5298 642*L 38544 0 5928 -642*L; -642*L -73*L^2 0 344*L^2 642*L -73*L^2;...
    0 0 5928 642*L 19272 -1458*L; 0 0 -642*L -74*L^2 -1458*L 172*L^2] ...
    %+ [zeros(4,6); zeros(2,4) [Mt St; St It]
    
%global stiffness equation
K2 = cK2*[24 6*L -24 6*L 0 0; 6*L 2*L^2 -6*L L^2 0 0;...
    -24 -6*L 48 0 -24 6*L; 6*L L^2 0 4*L^2 -6*L L^2;...
    0 0 -24 -6*L 24 -6*L; 0 0 6*L L^2 -6*L 2*L^2]

%reduced mass/stiffness equation
K4 = K2(3:6, 3:6)
M4 = M2(3:6, 3:6)


%pull out eigen values

[v, d] = eig(K4,M4)

%take square roots of diagonals

sorted  = sort(d)

%eig values in ascending order:
eigenvalues = sort(sorted(end, :))

freq = eigenvalues/(2*pi)

v1 = v(:, 4)
v2 = v(:,3)
v3 = v(:,2)
v4 = v(:,1)

v = [v1 v2 v3 v4]


%  ploteigenvector (L,ev,ne,nsub,scale);
% declare local variables here if required by language nv=ne*nsub+1; 
% Le=L/ne; dx=Le/nsub; k=1; x=v=zeroarray(nv); // declare and set to zero plot arrays 
% for (e=1,e<=ne,e++) // loop over elements
% 
% xi=Le*(e-1); vi=ev(2*e-1); qi=ev(2*e); vj=ev(2*e+1); qj=ev(2*e+2); 
% 
% for (n=1,n<=nsub,n++) %loop over subdivisions
% 
% xk=xi+dx*n; x=float(2*n-nsub)/nsub; // isoP coordinate vk=scale*(0.125*(4*(vi+vj)+2*(vi-vj)*(x^2-3)*x+
% 
% Le*(x^2-1)*(qj-qi+(qi+qj)*x))) // Hermitian interpolant k = k+1; x(k)=xk; v(k)=vk; // build plot functions
% 
% end % end n loop
% 
% end % end e loop
% 
% // plot v (vertical) vs x (horizontal) -- language dependent endprocedure











