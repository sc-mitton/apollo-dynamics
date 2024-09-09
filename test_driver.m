clear;
clc;

% Time (s)
dt = 0.01;      % Step size
tf = 30;        % Final time
t = 0:dt:tf;    % Time

% Initial conditions (deg and deg/s)
wx0 = 0;
wy0 = 0;
wz0 = 0;
psi0 = 0;
theta0 = 0;
phi0 = 0;

% Torques (N-m). You can change these torque vectors to validate your
% model and try out the different cases required in the project. When
% you send me your function, I will try out some torques to see if
% your model accurately predicts the response.
Mx = 176*cos(0.2*t);
%Mx = 0*ones(size(t))*0;

My = 54*ones(size(t));
%My = -0.9685 * ones(size(t))*0;

Mz = 98*sin(0.3*t);
%Mz = -0.4684 * ones(size(t))*0;

% Call the project function. Note that the initial values are passed in
% units of degrees and degrees/s, and the function returns solution
% vectors in units of degrees and degrees/s.
%
% ** CHANGE THIS TO THE NAME OF YOUR OWN FUNCTION. THE FUNCTION
% PARAMETERS AND ORDER SHOULD STAY THE SAME. THE NAME OF YOUR FUNCTION
% SHOULD BE lastname() ** Your function must return
% values of wx, wy, wz, psi, theta, and phi at each time step defined
% in the vector t.
[wx,wy,wz,psi,theta,phi]=mitton(wx0,wy0,wz0,psi0,theta0,phi0,t,Mx,My,Mz);

% disp('max wx:')
% disp(max(wx))
% disp('min wx:')
% disp(min(wx))
% disp('max wy:')
% disp(max(wy))
% disp('min wy:')
% disp(min(wy))
% disp('max wz:')
% disp(max(wz))
% disp('min wz:')
% disp(min(wz))
% disp('max psi:')
% disp(max(psi))
% disp('min psi:')
% disp(min(psi))
% disp('max theta:')
% disp(max(theta))
% disp('min theta:')
% disp(min(theta))
% disp('max phi:')
% disp(max(phi))
% disp('min phi:')
% disp(min(phi))

% Plot results
subplot(2,1,1);
plot(t,wx,t,wy,t,wz);
xlabel('t (s)');
ylabel('\omega (deg/s)');
%ylim([-.2,1.2])
legend('\omega_x','\omega_y','\omega_z','Location','northeastoutside');
set(gca,'FontSize',20)
subplot(2,1,2);
plot(t,psi,t,theta,t,phi);
%ylim([-40,1000])
xlabel('t (s)');
ylabel('\psi, \theta, \phi (deg)');
legend('\psi','\theta','\phi','Location','northeastoutside');
set(gca,'FontSize',20)