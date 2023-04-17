function PlotThinAirfoil(c,alpha,V_inf,p_inf,rho_inf,N)
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

%% Vortex Sheet (Gamma)
dx = c/N;
dx_center = dx/2;
for i = 1:N
    % distance to midpoint of panel
    z = dx*i - dx_center;
    % circulation at that location
    gamma = 2 * alpha * V_inf * sqrt( (1-z/c)/(z/c) );
    % circulation magnitude
    Gamma = gamma * dx;
    
    % stream line vortex around airfoil
    Psi_Gamma = (Gamma/(2*pi))*log(sqrt((x-z).^2 + y.^2));
    StreamFunction = StreamFunction + Psi_Gamma;
    
    % velocity potential vortex around airfoil
    Phi_Gamma = -(Gamma/(2*pi))*atan(y./(x-z));
    VelocityPotential = VelocityPotential + Phi_Gamma;
    
    % air velocities around airfoil
    Vx = Vx + ((Gamma*y)/(2*pi))./((x-z).^2+(y).^2);
    Vy = Vy - ((Gamma*(x-z))/(2*pi))./((x-z).^2+(y).^2);
end

%% Pressure Distribution
v = sqrt(Vx.^2+Vy.^2);
cP = 1-(v./V_inf).^2;
dP = 0.5 * rho_inf * V_inf^2;
P = cP * dP + p_inf;

%% Plot Stream Lines Psi
levmin = min(min(StreamFunction)); % defines the color levels 
levmax = max(max(StreamFunction));
levels = linspace(levmin,levmax,50)';

airfoilx = linspace(0,c,N);
airfoily = zeros(1,N);

figure(1)
subplot(2,2,1)
hold on
plot(airfoilx,airfoily,'k')
contour(x,y,StreamFunction,levels)
ylabel('y [m]')
xlabel('x [m]')
title('Stream Lines Thin Airfoil')
hold off

%% Plot Velocity Potential Phi
levmin = min(min(VelocityPotential)); % defines the color levels 
levmax = max(max(VelocityPotential));
levels = linspace(levmin,levmax,50)';

subplot(2,2,2)
hold on
plot(airfoilx,airfoily,'k')
contour(x,y,VelocityPotential,levels)
ylabel('y [m]')
xlabel('x [m]')
title('Velocity Potential Thin Airfoil')
hold off

%% Plot Pressure Contour Cp
levmin = min(min(P)); % defines the color levels 
levmax = max(max(P));
levels = linspace(levmin,levmax,50)';

subplot(2,2,3)
hold on
plot(airfoilx,airfoily,'k')
contour(x,y,P,levels)
ylabel('y [m]')
xlabel('x [m]')
title('Pressure Contours Thin Airfoil')
hold off

end