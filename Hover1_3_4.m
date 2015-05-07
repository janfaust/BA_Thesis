%Free Wake Method - Hover Version 1.3.4

 % - Aim: Complete for 1 Blade 

 % - Comment: Everything influences everything 


%initialization
close all;

%input variables for discretisation 
R=10;
Nr=16;
Nb=1; %Note: Functions only with one Nb atm, TODO: make adjustments to others 
t_rev=2; %time it takes for 1 revolution
omega=-(2*pi/t_rev);
t_final=1.5;
dt=0.2;
t_iter=ceil(t_final/dt)+1;

%wake input 
p_wake=zeros(3,t_iter,Nr); %(x,y,z), progresses with time, discretisation Nr 
Gamma_blade=zeros(Nr-1,1); %only one Blade for now
Gamma_wake=zeros(t_iter,Nr); %Gamma variable for Nr with time 
V=zeros(3,t_iter,Nr); %Velocity components for each point p_wake 
rc=0.5;

%rotational Matrix 
Rot_M=[cos(omega*dt) sin(omega*dt) 0; ...
        -sin(omega*dt) cos(omega*dt) 0; 0 0 1];
    
    r=linspace(0.1*R,R,Nr); %seperation along r 
    p=zeros(3,Nr,Nb); %p_blade decleration 

%blade configuration, initialization 
    for k=1:Nb
        p(1,:,k)=r*cos((2*pi/Nb)*(k-1));
        p(2,:,k)=r*sin((2*pi/Nb)*(k-1));
    end
%Gamma_blade (simplification: Gamma=-4*pi*sin(alpha)*q_inf
sin_alpha=sin(10*pi/180); %assumption 10 deg constant along r, TODO: subprogram variable alpha(r)
q_inf=zeros(Nr,1); 
q_inf=r*abs(omega); %flow velocity component (only tangential to blade atm)
Gamma_temp=-4*pi*sin_alpha*q_inf;
for ii=1:Nr-1
Gamma_blade(ii)=0.5*(Gamma_temp(ii)+Gamma_temp(ii+1)); %constant gamma distribution along Blade 
end

tic;
% Simulation
for ii=1:1:t_iter
   
        %shed wake  
        for tt=2:Nr-1  %for each timestep the wake is shed, based on previous gamma_blade along r 
            Gamma_wake(ii,1)=-Gamma_blade(1); 
            Gamma_wake(ii,tt)=Gamma_blade(tt-1)-Gamma_blade(tt); 
            Gamma_wake(ii,Nr)=Gamma_blade(Nr-1);
        end
        
        for jj=1:Nr         
    %position, coordinates 
    p(:,jj,k)=Rot_M*p(:,jj,k); %new blade position
    p_wake(:,ii,jj)=p(:,jj,k); %new wake position 
   
        end
    % ----------------------------------------------------------------%
    %BiotSavart - Euler 
    %function [ V ] = BiotSavartLaw3( ra,rb,rp,Gamma )
    %p_wake=zeros(3,t_iter,Nr);
    %Gamma_wake=zeros(t_iter,Nr); 
    %V=zeros(3,t_iter,Nr);
    if (ii>1) %only calculate when timestep is bigger than 1 
    V(:,:,:)=0;   %clear V before each new iteration 
     for kk=1:ii-1 
         
   %tip vortex -> other vorticies 
           for nn=1:Nr-1
                 V(:,kk,nn)=sum(BiotSavartLaw3(p_wake(:,kk,Nr),p_wake(:,kk+1,Nr),...
                 p_wake(:,kk,nn),Gamma_wake(kk,end)));
           end 
   %tip vortex ->self inducing effect only 
             if kk>1
                for ll=1:kk-1    
                V(:,ll,Nr)=sum(BiotSavartLaw3(p_wake(:,kk,Nr),p_wake(:,kk+1,Nr),...
                p_wake(:,ll,Nr),Gamma_wake(kk,end)));
                    if ll>2
                        for mm=1:ll-2
                            V(:,ll,Nr)=sum(BiotSavartLaw3(p_wake(:,mm,Nr),p_wake(:,mm+1,Nr),...
                            p_wake(:,ll,Nr),Gamma_wake(kk,end)));
                        end 
                    end
                end   
             end
        
     end
   % ----------------------------------------------------------------%
    %PROBLEM AREA: Velocity self inducing tip vortex
    
   %  for kk=1:ii-1 
   %  V(:,1,Nr)=sum(BiotSavartLaw3(p_wake(:,kk,Nr),p_wake(:,kk+1,Nr),...
   %    p_wake(:,1,Nr),Gamma_wake(kk,end)));
   % end
    %problem wenn p_wake and p_gamma are on the same point -> NaN
    %reason being: 
    % ----------------------------------------------------------------%
    
    %Euler Equation 
    %p_wake()=p_wake()+dt*V(); 
    %Note: vi is only used for displaying purposes, from Momentum Theory
    %assumption: vi is added to every component, as a constant 
    %in hover, vi=sqrt(Thrust/(2*rho*pi*R^2)), Thrust=m*g
   
    vi=zeros(3,1);
    vi(3)=0.1;
    
    for oo=1:Nr
        for pp=1:ii-1
            p_wake(:,pp,oo)=p_wake(:,pp,oo)+dt*V(:,pp,oo)-vi(:);
        end
    end
    
    % ----------------------------------------------------------------%
    end %if ii>1
    %counter:iteration
       iteration=ii
       
end
toc;
%plot end result 
figure(1)
%plot near wake using plot3
for ii=1:Nr
plot3(p_wake(1,:,ii),p_wake(2,:,ii),p_wake(3,:,ii));
hold on; 
end

%plot current blade position
plot3(p(1,:,:),p(2,:,:),p(3,:,:),'k','LineWidth',3)
axis([-2*R 2*R -2*R R*2 -10 0]);
ylabel('y');
xlabel('x');
zlabel('z');
title('Wake for hover flight, one blade')

