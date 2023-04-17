function CL = PlotThickAirfoilCl(c,alpha,V_inf,p_inf,rho_inf,N,xx)
%% Plotting domain/range
DR = [-1 3 -1 1]; %[xmin xmax;ymin ymax] for plotting axes

% Grid of points over domain/range
[x,y] = meshgrid(linspace(DR(1),DR(2),100),linspace(DR(3),DR(4),100));

%% Flow Parameters
% uniform stream function
Psi_inf = V_inf * (y*cos(alpha) - x*sin(alpha));
StreamFunction = Psi_inf;
% uniform velocity potential
Phi_inf = V_inf * (x*cos(alpha) - y*sin(alpha));
VelocityPotential = Phi_inf;
% uniform velocity
Vx = V_inf * cos(alpha)*ones(100);
Vy = V_inf * sin(alpha)*ones(100);

%% Airfoil Shape
t = xx/100;
dx = c/N;
XB = zeros(1,N);
YB = zeros(1,N);
for i = 1:(N/2)
    % distance to midpoint of panel
    z = 2*dx*(i-1);
    % thickness of airfoil from chord
    yt = (t*c/0.2)*(0.2969*sqrt(z/c)-0.1260*(z/c)-0.3516*(z/c)^2+0.2843*(z/c)^3-0.1036*(z/c)^4);
    % XB is vector of airfoil X points
    XB(i) = z;
    XB(N-i+1) = z;
    % YB is vector of corresponding Y points
    YB(i) = yt;
    YB(N-i+1) = -1*yt;
end
% figure(2)
% plot(XB,YB)

%% Call Vortex Panel
% X,Y are midpoints of panels
[CL,CP,GAMMA,X,Y] = Vortex_Panel(XB,YB,V_inf,alpha,0);
CL = (pi*alpha*rho_inf*c*V_inf^2)/()
end