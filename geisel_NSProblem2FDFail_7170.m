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
%--------------------------------------------------------------------------

%Establish parameters

ymesh=125; %number of evenly spaced y values to be used in the approximation
xmesh=65; %number of evenly spaced y values to be used in the approximation
w=2; %half of the width of the surface, perpendicular to the flow.
inty=2*w ;%interval length in the y direction
leng=4*w; %length of the ramp in the x direction
h=inty/(ymesh-1); %diffence between each value of y
eta=0.005 ;%1/Reynold's number /Assume Reynolds number of 200 in this case. 
P=-1 ;%Constant for change in Pressure in the positive x direction (deltaP/deltax)
        %P must be negative since the the water flows to lower pressure.
rho=0.9983; %density of water at 20 C       
y=linspace(-w,w,ymesh);
theta=20; %angle of declination measured in degrees
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
M=zeros(ymesh,rmesh);
Mold=zeros(ymesh,1);
Mold(2:ymesh-1,1)=1;
const=-P+rho*gr; %repeated constants in calculation
Mnew=zeros (ymesh,1);
Mnew(2:ymesh-1,1)=1;
for i=1:50
    Mhold=Mnew;
    for j=2:ymesh-1
    leadcoef=1/(((rho*Mhold(j))/(2*k))-(eta/(k^2)));
    midcoef=(((rho*Mhold(j))/(2*k))+(eta/(k^2)));
    Mnew(j)=leadcoef*(const+(eta/(k^2))*Mhold(j-1)+midcoef*Mold(j)-((4*eta)/(k^2))*Mhold(j)+(eta/(k^2))*Mhold(j+1));
    end
    M(:,i)=Mhold;
    Mold=Mhold;
end

surf(M)