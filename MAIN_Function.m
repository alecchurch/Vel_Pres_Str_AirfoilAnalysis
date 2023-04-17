% Alec Church
% CA 2 Flow over Airfoils
% 2/27/2022
% The purpose of this assignment is to study the effects the shape of an
% airfoil has on the stream lines, velocity potentials, and coefficients of
% pressure on the surrounding air.

%% Problem 1: Visualizing flow around thin airfoils
c = 2;                  %Chord length (m)
ALPHA = 9 * (pi/180);   %Angle of attack (radians)
VINF = 60;              %Freestream velocity (m/s)
p_inf = 85.5 * 10^3;    %Freestream pressure (Pa)
rho_inf = 1;            %Air density (kg/m^3)

N = 100;                %Number of panels
%% Problem 1 Function PlotThinAirfoil.m
PlotThinAirfoil(c,ALPHA,VINF,p_inf,rho_inf,N);
fprintf('Manually varying the number of panels, N, between 10 and 10,000 showed \nnoticable accuracy at 100 panels and beyond. Noncontinuous steam \nlines, velocity potentials, and pressure contours were present until \nthe number of panels reached N = 100. At N > 100 there was no noticable \nincrease in performance or accuray. \n')

%% Problem 2 Function PlotThickAirfoil.m
% NACA 0012 symmetric airfoil
N=500;
xx = 12;
PlotThickAirfoil(c,ALPHA,VINF,p_inf,rho_inf,N,xx)
fprintf('Manually varying the number of panels, N, up to 1000 showed\ngreat accuracy at 1000 panels. After passing that number, the next \norder of magnitude, 10,000 panels, was inefficient to run and had no \ngreater increase in accuracy. \n ')
