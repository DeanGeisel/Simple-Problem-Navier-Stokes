%--------------------------------------------------------------------------
%For: 'A Brief Exploration of the Navier Stokes Equations'
%By: Dean Geisel
%MTH7170 Final Project
%Instructor: Dr. Sulman
%--------------------------------------------------------------------------
%The purpose of this code is to approximate the velocity field of a fluid
%crossing a flat surface, with zero velocity side walls parallel to the
%flow. The assumptions are: No slip sides, Horizontal surface, fully
%developed flow, steady state, incompressible fluid
%There is velocity only in the x direction, so only u is solved for.
%--------------------------------------------------------------------------
m=[9 33 65 129];
%for i=1:4
%establish parameters
mesh=m(3); %number of y values to be 
w=1; %half of the width of the surface, perpendicular to the flow.
inty=2*w ;%interval length in the y direction
h=inty/(mesh-1); %diffence between each value of y
eta=0.005 ;%1/Reynold's number /Assume Reynolds number of 200 in this case.
P=-0.5 ;%Constant for change in Pressure in the positive x direction (deltaP/deltax)
        %P must be negative since the the water flows to lower pressure.
y=linspace(-w,w,mesh);
x=linspace(0,4*w,2*mesh);
%--------------------------------------------------------------------------

%create coefficeint matrix a for the approximation

A=zeros(mesh-2); %establish empty A ready for diagonals
e=ones(mesh-2,1); %vector for super and sub diagonals
d=-2*e; %Vector as negative double e for main diagonal
A=spdiags([e d e],[-1 0 1],mesh-2,mesh-2); %fill A with the diagonal entries
A=1/(h^2)*A;
%--------------------------------------------------------------------------

%Create Vector F "right side"
%No subtractions because of the no slip conditions giving first and last
%values as zero.

F=(P/eta)*e;

%--------------------------------------------------------------------------

%Solve for velocity profile

Ushort=A\F;

U=zeros(mesh,1);
U(2:mesh-1)=Ushort;

%--------------------------------------------------------------------------

%Extend for diagram

Ucont=zeros(mesh,2*mesh);
Ucont=repmat(U,[1,2*mesh]);

%--------------------------------------------------------------------------

%Plot surfaces

figure(1)
plot (U,y)
figure(2)
surf(Ucont);

%end
