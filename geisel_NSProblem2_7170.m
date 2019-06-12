%--------------------------------------------------------------------------
%For: 'A Brief Exploration of the Navier Stokes Equations'
%By: Dean Geisel
%MTH7170 Final Project
%Instructor: Dr. Sulman
%--------------------------------------------------------------------------
%The purpose of this code is to approximate the velocity field of a fluid
%crossing a tilted surface at angle theta, with zero velocity side walls parallel to the
%flow. The assumptions are: Laminar Flow, No slip sides, angle of declination theta, non
%fully developed flow,steady state, incompressible fluid.
%This code follows the method of rotating space so that the flow is
%horizontal, excluding calculations for the z direction.

%The method used is Newton's Method for nonlinear PDE's

%This code produced incorrect data 

%--------------------------------------------------------------------------

%Establish parameters

ymesh=5; %number of evenly spaced y values to be used in the approximation
xmesh=9; %number of evenly spaced x values to be used in the 
N=ymesh*xmesh; %abbreviation for common used number
w=4; %half of the width of the surface, perpendicular to the flow.
inty=2*w ;%interval length in the y direction
leng=4*w ;%length of the ramp in the x direction
h=inty/(ymesh-1); %diffence between each value of y also used for r
eta=0.005 ;%1/Reynold's number /Assume Reynolds number of 200 in this case. 
P=0 ;%Constant for change in Pressure in the positive x direction (deltaP/deltax)
        %P must be negative since the the water flows to lower pressure.
rho=0.9983 ;%density of water at 20 C       
y=linspace(-w,w,ymesh);
theta=10; %angle of declination measured in degrees
g=9.8; %force of gravity in m/second^2

%--------------------------------------------------------------------------

%calculations to rotate figure
    
    %gravity
gr=(sind(theta))*g;  %gravity in r direction (along the ramp, effective gravity)
gz=(sind(theta))*gr; %effective gravity in z direction
gx=(cosd(theta))*gr; %effective gravity in x direction

rmesh=xmesh; %keep the number of intervals the same, even if length is longer
k=h; %delta ramp is the same as the change in y to facilitate calculation

%--------------------------------------------------------------------------


U=zeros(N,1);

for j=1:4

UM=reshape(U,ymesh,xmesh);


%forming F-----------------------------------------------------------------

%F body

FM=zeros(ymesh,xmesh); %Matrix version of F to match UM
C=(-P+rho*gr)*k^2; %constant independent of location
for p=2:xmesh-1
    for q=2:ymesh-1
        L=rho*k*UM(q,p)*(UM(q,p+1)-UM(q,p));
        R=eta*(UM(q,p-1)+UM(q-1,p)-4*UM(q,p)+UM(q+1,p)+UM(q,p+1));
      FM(q,p)=L-C-R ;
    end
end

%F left border (initial x, not corners)
 
 for q=2:ymesh-1
    Uglft=2*UM(q,1)-UM(q,2);
    L=rho*k*UM(q,1)*(UM(q,2)-UM(q,1));
    R=eta*(Uglft+UM(q-1,1)-4*UM(q,1)+UM(q+1,1)+UM(q,2));
    FM(q,1)=L-C-R;
 end
 
 %F upper border (not corners at 1,1)
 
 for p=2:xmesh-1
    Ugupp=2*UM(1,p)-UM(2,p);
    L=rho*k*UM(1,p)*(UM(1,p+1)-UM(1,p));
      R=eta*(UM(1,p-1)+Ugupp-4*UM(1,p)+UM(2,p)+UM(1,p+1));
    FM(1,p)=L-C-R;
 end
 
 %F lower border (not corners at 1,1)
 
 for p=2:xmesh-1
    Uglow=2*UM(ymesh,p)-UM(ymesh-1,p);
    L=rho*k*UM(ymesh,p)*(UM(ymesh,p+1)-UM(ymesh,p));
    R=eta*(UM(ymesh,p-1)+UM(ymesh-1,p)-4*UM(ymesh,p)+Uglow+UM(ymesh,p+1));
    FM(ymesh,p)=L-C-R;
 end
 
 %F right border (final x, not corners)
 
 for q=2:ymesh-1
    Ugrht=2*UM(q,xmesh)-UM(q,xmesh-1);
    L=rho*k*UM(q,xmesh)*(Ugrht-UM(q,xmesh));
    R=eta*(UM(q,xmesh-1)+UM(q-1,xmesh)-4*UM(q,xmesh)+UM(q+1,xmesh)+Ugrht);
    FM(q,xmesh)=L-C-R;
 end
 
 %Corners of F
 %Upper left
 L=rho*k*UM(1,1)*(UM(1,2)-UM(1,1));
 R=eta*(2*UM(1,1)-UM(1,2)+2*UM(1,1)-UM(2,1)-4*UM(1,1)+UM(2,1)+UM(1,2));
 FM(1,1)=L-C-R;
 
 %Lower Left
 L=rho*k*UM(ymesh,1)*(UM(ymesh,2)-UM(ymesh,1));
 R=eta*(2*UM(ymesh,1)-UM(ymesh,2)+UM(ymesh-1,1)-4*UM(ymesh,1)+2*UM(ymesh,1)-UM(ymesh-1,1)+UM(ymesh,2));
 FM(ymesh,1)=L-C-R;
 
 %Upper Right
  L=rho*k*UM(1,xmesh)*(2*UM(1,xmesh)-UM(1,xmesh-1)-UM(1,xmesh));
  R=eta*(UM(1,xmesh-1)+2*UM(1,xmesh)-UM(2,xmesh)-4*UM(1,xmesh)+UM(2,xmesh)+2*UM(1,xmesh)-UM(1,xmesh-1));
      FM(1,xmesh)=L-C-R ;
 
  %Lower Right
  L=rho*k*UM(ymesh,xmesh)*(2*UM(ymesh,xmesh)-UM(ymesh,xmesh-1)-UM(ymesh,xmesh));
   R=eta*(UM(ymesh,xmesh-1)+UM(ymesh-1,xmesh)-4*UM(ymesh,xmesh)+2*UM(ymesh,xmesh)-UM(ymesh-1,xmesh)+2*UM(ymesh,xmesh)-UM(ymesh,xmesh-1));
      FM(ymesh,xmesh)=L-C-R ; 
  
 
 F=reshape(FM,N,1);
 
%Forming the Jacobian Main diagonal

Jhold=zeros(ymesh,xmesh);
for p=1:xmesh-1
    for q=1:ymesh
        Jhold(q,p)=rho*k*(UM(q,p+1)-2*UM(q,p))+4*eta;
    end
end

for q=1:ymesh
    Ugst=2*UM(q,xmesh)-UM(q,xmesh-1);
    Jhold(q,xmesh)=rho*k*(Ugst-2*UM(q,p))+4*eta;
end


JV=reshape(Jhold,N,1);

%Forming the Right diagonal dependent upon U
Rhold=zeros(ymesh,xmesh-1);
for p=1:xmesh-1
    for q=1:ymesh
        Rhold(q,p)=rho*k*(UM(q,p))-eta;
    end
end

RV=reshape(Rhold,ymesh*(xmesh-1),1);

%Assembling the Jacobian

I=eye(xmesh);
IJV=spdiags(JV,0,N,N);
ILV=spdiags(RV,-ymesh,N,N);
e=ones(ymesh,1);
ek=ones(xmesh,1);
K=spdiags(ek,-1,xmesh,xmesh);
T=spdiags([-eta*e -eta*e],[-1 1],ymesh,ymesh);
Er=spdiags(-eta*e,0,ymesh,ymesh);
J=(kron(I,T)+kron(K,Er)+IJV+ILV');
%figure (1)
%spy(J)

%Solve for Delta U
DU=J\-F;

%Solving for U(n+1)
Un=DU+U;
UnM=reshape(Un,ymesh,xmesh);
UnM(:,1)=0;
UnM(1,:)=0;
UnM(ymesh,:)=0;
Un=reshape(UnM,N,1);
U=Un;

end

%Test=reshape(Inthold(:,2),ymesh,xmesh)



UM=reshape(U,ymesh,xmesh);  
  
figure(2)
surf(UM)
        