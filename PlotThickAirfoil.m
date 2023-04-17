function PlotThickAirfoil(c,alpha,V_inf,p_inf,rho_inf,N,xx)
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gamma(i) at  X(i) Y(i)

%% Vortex Sheet (Gamma)
for i = 1:length(GAMMA)
    x_mid = X(i);
    y_mid = Y(i);
    Gamma = GAMMA(i);
    
    % stream line vortex around airfoil
    Psi_Gamma = (Gamma/(2*pi)).*log(sqrt((x-x_mid).^2 + (y-y_mid).^2));
    StreamFunction = StreamFunction + Psi_Gamma;
    
    % velocity potential vortex around airfoil
    Phi_Gamma = -(Gamma/(2*pi)).*atan((y-y_mid)./(x-x_mid));
    VelocityPotential = VelocityPotential + Phi_Gamma;
    
    % air velocities around airfoil
    Vx = Vx + ((Gamma*(y-y_mid))/(2*pi))./((x-x_mid).^2+(y-y_mid).^2);
    Vy = Vy - ((Gamma*(x-x_mid))/(2*pi))./((x-x_mid).^2+(y-y_mid).^2);
end

%% Pressure Distribution
v = sqrt(Vx.^2+Vy.^2);
cP = 1-(v./V_inf).^2;
dP = 0.5 * rho_inf * V_inf^2;
P = cP .* dP + p_inf;

%% Plot Stream Lines Psi
levmin = min(min(StreamFunction)); % defines the color levels 
levmax = max(max(StreamFunction));
levels = linspace(levmin,levmax,50)';

figure(2)
subplot(2,2,1)
hold on
plot(XB,YB,'k')
contour(x,y,StreamFunction,levels)
ylabel('y [m]')
xlabel('x [m]')
title('Stream Lines NACA 0012 Airfoil')
hold off

%% Plot Velocity Potential Phi
levmin = min(min(VelocityPotential)); % defines the color levels 
levmax = max(max(VelocityPotential));
levels = linspace(levmin,levmax,50)';

subplot(2,2,2)
hold on
plot(XB,YB,'k')
contour(x,y,VelocityPotential,levels)
ylabel('y [m]')
xlabel('x [m]')
title('Velocity Potential  NACA 0012 Airfoil')
hold off

%% Plot Pressure Contour Cp
levmin = min(min(P)); % defines the color levels 
levmax = max(max(P));
levels = linspace(levmin,levmax,50)';

subplot(2,2,3)
hold on
plot(XB,YB,'k')
contour(x,y,P,levels)
ylabel('y [m]')
xlabel('x [m]')
title('Pressure Contours  NACA 0012 Airfoil')
hold off

end