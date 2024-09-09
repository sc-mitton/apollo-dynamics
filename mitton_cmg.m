function [wx,wy,wz,psi,theta,phi] = mitton(wx0,wy0,wz0,psi0,theta0,phi0,t,th1,th2,th3)
    
    % Inertial Properties %%
    %%%%%%%%%%%%%%%%%%%%%%%%

    % Command Module 
    weight_cm = 9730/2.20462442018; % kg
    COM_cm = [1043.1;-.1;7.8] ./39.3701; % center of mass in A frame in m
    I_cm = [4474,0,0;0,3919,0;0,0,3684] * 1.3558179619; % I in kg*m^2 about the COM of each component

    % Service Module
    weight_sm = 9690/2.20462442018;
    COM_sm = [908.2;.7;-.6] ./39.3701; 
    I_sm= [6222,0,0;0,10321,0;0,0,10136] * 1.3558179619;

    % Propellant
    weight_prop = 37295/2.20462442018;
    COM_prop = [905.9;5.6;-2.4] ./39.3701;
    I_prop = [19162,0,0;0,19872,0;0,0,26398]*1.3558179619;

    % Command Service Module (including propellants)
    weight_csm = weight_cm+weight_sm+weight_prop;

    % about A frame
    for i = 1:3
        COM_csm(i,1) = (COM_sm(i)*weight_sm + COM_cm(i)*weight_cm + COM_prop(i)*weight_prop ) / (weight_csm);
    end

    d = -COM_sm + COM_csm;
    I_sm_B  = I_sm + weight_sm*(d'*d*eye(3) - d*d');

    d = -COM_cm + COM_csm;
    I_cm_B = I_cm + weight_cm*(d'*d*eye(3) - d*d');

    d = -COM_prop + COM_csm;
    I_prop_B = I_prop + weight_prop*(d'*d*eye(3) - d*d');

    % total Inertia matrix including propellant
    I_csm_B = (I_sm_B + I_cm_B + I_prop_B);
    
    COM_sm_simp=COM_sm;
    COM_sm_simp(2:3) = 0;
    COM_cm_simp = COM_cm;
    COM_cm_simp(2:3) = 0;
    COM_prop_simp = COM_prop;
    COM_prop_simp(2:3) = 0;

    % about A frame
    for i = 1:3
        COM_csm_simp(i,1) = (COM_sm_simp(i)*weight_sm + COM_cm_simp(i)*weight_cm + COM_prop_simp(i)*weight_prop ) / (weight_csm);
    end

    d = -COM_sm_simp + COM_csm_simp;
    I_sm_B_simp  = I_sm + weight_sm*(d'*d*eye(3) - d*d');

    d = -COM_cm_simp + COM_csm_simp;
    I_cm_B_simp = I_cm + weight_cm*(d'*d*eye(3) - d*d');

    d = -COM_prop_simp + COM_csm_simp;
    I_prop_B_simp = I_prop + weight_prop*(d'*d*eye(3) - d*d');

    % total Inertia matrix including propellant
    I_csm_B_simp = I_sm_B_simp + I_cm_B_simp + I_prop_B_simp;
        
    % Task D.b
    %%%%%%%%%%%%%%%%%%%%%%
    
    % Find: minimum rotor inertia Ir needed to achieve desired angle
    % accelerations of wx_dot = 30 deg/s^2 , wy_dot = wz_dot = 2.5 deg/s^2

%     Omega = 6600; % rpm
%     Omega = 6600 *(1/60) * 2*pi ; % rad/s
%     Omega1 = Omega; Omega2 = Omega; Omega3 = Omega;
%     
%     thetad_max = 30; % deg/s
%     thetad_max = 30*pi/180; % rad/s
%     thd1 = thetad_max; thd2 = thetad_max; thd3 = thetad_max;
% 
%     wxd = -5*pi/180*0;
%     wyd = -2.5*pi/180*0;
%     wzd = -2.5*pi/180;
%     omd = [wxd;wyd;wzd];
%     th1 = 0;
%     th2 = 0;
%     th3 = 0;
%     
%     I = I_csm_B_simp;
    
%     M = [-cos(th3)*thd3*Ir*Omega3 + sin(th1)*thd1*Ir*Omega1;
%          -cos(th1)*thd1*Ir*Omega1 + sin(th2)*thd2*Ir*Omega2;
%          -cos(th2)*thd2*Ir*Omega2 + sin(th3)*thd3*Ir*Omega3;];
%     
%     I(1,1)*omd(1) - I(1,2)*omd(2) - I(1,3)*omd(3) = M(1);
%     I(2,2)*omd(2) - I(2,3)*omd(3) - I(1,2)*omd(1) = M(2);
%     I(3,3)*omd(3) - I(1,3)*omd(1) - I(2,3)*omd(2) = M(3);
    
%     Ir = (I(1,1)*omd(1) - I(1,2)*omd(2) - I(1,3)*omd(3)) / (-cos(th3)*thd3*Omega3 + sin(th1)*thd1*Omega1);
%     disp(Ir)
%     Ir = (I(2,2)*omd(2) - I(2,3)*omd(3) - I(1,2)*omd(1)) / (-cos(th1)*thd1*Omega1 + sin(th2)*thd2*Omega2);
%     disp(Ir)
%     Ir = (I(3,3)*omd(3) - I(1,3)*omd(1) - I(2,3)*omd(2)) / (-cos(th2)*thd2*Omega2 + sin(th3)*thd3*Omega3);
%     disp(Ir)
    
    % Task E.a
    %%%%%%%%%%%%%%%%%%%%%%
    
    % I found that the minimum Ir needed is 11.8928 kg*m^2
    Ir = 11.8928;
    Omega = 6600*(1/60)*2*pi; % rad/s
    Omega1 = Omega; Omega2 = Omega; Omega3 = Omega;
    
    % calculate derivatives for thetas
    th1 = th1 * pi/180;
    th2 = th2 * pi/180;
    th3 = th3 * pi/180;
    
    dt = t(2) - t(1);
    thd1 = diff(th1)/dt;
    thd2 = diff(th2)/dt;
    thd3 = diff(th3)/dt;
    
    thd1(end+1) = thd1(end);
    thd2(end+1) = thd2(end);
    thd3(end+1) = thd3(end);
    
    % moments enacted on the ship
    Mx = -cos(th3).*thd3*Ir*Omega3 + sin(th1).*thd1*Ir*Omega1;
    My = -cos(th1).*thd1*Ir*Omega1 + sin(th2).*thd2*Ir*Omega2;
    Mz = -cos(th2).*thd2*Ir*Omega2 + sin(th3).*thd3*Ir*Omega3;
    
    initials = [wx0,wy0,wz0,psi0,theta0,phi0];
    initials = initials*pi/180;
    
    times = t;
    [t,y1] = ode45(@(t,x) spencer(t,x,Mx,My,Mz,I_csm_B_simp,times),t,initials);
    
    y = y1*180/pi;
    
    wx = y(:,1);
    wy = y(:,2);
    wz = y(:,3);
    psi = y(:,4);
    theta = y(:,5);
    phi = y(:,6);

end

function xdot = spencer(t,x,Mx,My,Mz,I,times) 

    om1 = x(1);
    om2 = x(2);
    om3 = x(3);
    theta = x(5);
    phi = x(6);
    
    wx = om1;
    wy = om2;
    wz = om3;
    
    Mx = interp1(times,Mx,t);
    My = interp1(times,My,t);
    Mz = interp1(times,Mz,t);
    
    % EOM
    B = [I(1,2)*wx*wz - I(1,3)*wx*wy - (I(2,2)-I(3,3))*wy*wz - I(2,3)*(wy^2 - wz^2);
         I(2,3)*wx*wy - I(1,2)*wy*wz - (I(3,3)-I(1,1))*wx*wz - I(1,3)*(wz^2 -wx^2);
         I(1,3)*wy*wz - I(2,3)*wx*wz - (I(1,1)-I(2,2))*wx*wy - I(1,2)*(wx^2 - wy^2)];
    wdot = I\([Mx;My;Mz] - B);

    om1dot = wdot(1);
    om2dot = wdot(2);
    om3dot = wdot(3);
    
    xdot = [om1dot;
            om2dot;
            om3dot;
            (1/cos(theta))*(om2*sin(phi) + om3*cos(phi));
            om2*cos(phi) - om3*sin(phi);
            (1/cos(theta))*(om2*sin(theta)*sin(phi) + om3*sin(theta)*cos(phi)) + om1];
end