function [ V ] = BiotSavartLaw3( ra,rb,rp,Gamma)
%Induced velocity calculated using the BiotSavartLaw (ra,rb,rp,gamma,rc),
% Note: Added correction, p.601 equation 10.48 
% Function used for implimentation of BiotSavartLaw
% ra=starting point of vortex element 
% rb=ending point of vortex element
% rp=point velocity is being calculated for
% Gamma=vorticity strength of vortex element 
% rc=vortex core radius 
n=1.2;
rc=0.5;

r1=rp-ra;
r2=rp-rb;
l12=rb-ra;

cos_theta1=(l12'*r1)/(sqrt(sum(l12(:).^2))*sqrt(sum(r1(:).^2)));
cos_theta2=(l12'*r2)/(sqrt(sum(l12(:).^2))*sqrt(sum(r2(:).^2)));

h=sqrt(sum(r1(:).^2))*sin(acos(cos_theta1));
%h=r1*sin(theta1)=r2*sin(theta2)

kreuzprod=cross(l12,r1);
e=kreuzprod/(sqrt(sum(kreuzprod(:).^2)));

%V=(Gamma/(4*pi*h))*e;
V=(Gamma/(4*pi))*(h/(rc^(2*n)+h^(2*n))^(1/n))*(cos_theta1-cos_theta2)*e;

end

