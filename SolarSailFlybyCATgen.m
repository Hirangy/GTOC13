clear;
clc;
data = readmatrix('/Users/hirangy/Downloads/gtoc13_planets.csv', 'NumHeaderLines', 1);

%Filter Planet Ratios
aratios = [1,2,2];
for i = 1:6
    for j = 1:7
        ratio = data(j,5)/data(i,5);
        if ratio ~= 1 && ratio < 6 && ratio > 0.6 && i ~= 5
            aratios = [aratios; ratio,i,j];
        end
    end
end

%INITIALIZATION
doPLOT = 0;
mu = 139348062043.343; %[km^3/s^2]
r0 = 149597870.691; %[km]
Y = 365.25*86400;

tstart = 0*Y;
RK_steps = 200;

flyby_pla = 2;
target_pla = 3;
target = 3;
a = data(flyby_pla,5);
e = 0.001;
inc = 0*pi/180;
RAAN = 0*pi/180;
omega = 0*pi/180;
M0 = 0*pi/180;

T = 2*pi*sqrt(a^3/mu);
M = mod(M0 + 2*pi*tstart/T, 2*pi);
E = M;
f = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
r = a.*(1-e.^2)./(1+e.*cos(f));
v = sqrt(2*mu./r-mu/a);

RV_flyby(1:3,1) = [r.*(cos(f+omega)*cos(RAAN)-sin(f+omega)*cos(inc)*sin(RAAN)),...
    r.*(cos(f+omega)*sin(RAAN)+sin(f+omega)*cos(inc)*cos(RAAN)),...
    r.*(sin(f+omega)*sin(inc))];
gamma = atan(e.*sin(f)./(1+e.*cos(f)));
RV_flyby(4:6,1) = [v.*(-sin(f+omega-gamma)*cos(RAAN)-cos(f+omega-gamma)*cos(inc)*sin(RAAN)),...
    v.*(-sin(f+omega-gamma)*sin(RAAN)+cos(f+omega-gamma)*cos(inc)*cos(RAAN)),...
    v.*(cos(f+omega-gamma)*sin(inc))];

v_o = sqrt(mu/a);

v_in_vals = (0.05:0.005:0.40)*v_o;

n_DeltaM = 360;
DeltaM_vals = linspace(0,(2*pi)*(n_DeltaM-1)/n_DeltaM,n_DeltaM);

t_ahead_vals = 0.1:0.01:4;

n_gamma = 180;
gamma_vals = linspace(-pi,pi,n_gamma);

alpha_vals = [-0.5,-0.25,0,0.25,0.5,pi/2];
delta_vals = zeros(1,numel(alpha_vals));

n_coast = 1;
alpha_coast = [pi/2];
delta_coast = [0];

CATcell = cell(numel(v_in_vals), 1); % outer layer, one per iv
for iv = 1:numel(v_in_vals)
    CATcell{iv} = cell(n_DeltaM, size(aratios,1)); % inner 2D cell array
end

CAT = struct('vi_vals',v_in_vals/v_o,...
    'DeltaM_vals',DeltaM_vals,...
    'gamma_vals',gamma_vals,...
    'alpha_coast',alpha_coast(1:n_coast),...
    'delta_coast',delta_coast(1:n_coast),...
    'aratios',aratios,...
    'cell',CATcell...
    );

%Planned duration of sailing (fraction of Y)

%Original SC path
f_vals = (0:0.01*pi:2*pi)';

sail_max = 0.5;

for i = 1:5
    a_pla = data(i,5);
    e_pla = 0;
    inc_pla = 0*pi/180;
    RAAN_pla = 0*pi/180;
    omega_pla = 0*pi/180;
    r_path = a_pla.*(1-e_pla.^2)./(1+e_pla.*cos(f_vals));
    R_path_pla{i} = [r_path.*(cos(f_vals+omega_pla)*cos(RAAN_pla)-sin(f_vals+omega_pla)*cos(inc_pla)*sin(RAAN_pla)),...
        r_path.*(cos(f_vals+omega_pla)*sin(RAAN_pla)+sin(f_vals+omega_pla)*cos(inc_pla)*cos(RAAN_pla)),...
        r_path.*(sin(f_vals+omega_pla)*sin(inc_pla))];
end

plotrange1 = 2*r0;

tic;
%FIRST STAGE: sweep large arc to catch any planets
if doPLOT
    figure('Position', [100 100 1100 700]);
end

Gamma_vec = [sin(gamma_vals);cos(gamma_vals);zeros(1,numel(gamma_vals))];
delete(gcp('nocreate'));   % stop all workers
parpool(8);
parfor iv = 1:numel(v_in_vals)
    RV_out = zeros(6,n_gamma);
    v_in = v_in_vals(iv);
    for i = 1:n_gamma
        RV_out(:,i) = [RV_flyby(1:3,1);RV_flyby(4:6,1)+Gamma_vec(:,i)*v_in];
    end

    caught = zeros(1,n_DeltaM*size(aratios,1));

    COE_out = RV2COE(RV_out,mu);
    a_out = COE_out(1,:);
    T_out = 2*pi*sqrt(a_out.^3/mu);


    for k = 1:numel(t_ahead_vals)
        t_ahead = t_ahead_vals(k);
        tstep = t_ahead*T/RK_steps;
        t_vals = tstart:tstep:tstart+t_ahead*T;
        R_pla = PlaPos(tstart+t_ahead*T,data,0);
        Rfinal1 = zeros(3,numel(alpha_vals));

        R_pla_target = zeros(n_DeltaM*size(aratios,1),3);

        for im = 1:n_DeltaM
            for n = 1:size(aratios,1)
                F = 2*pi*t_ahead /aratios(n,1)^1.5;
                R_pla_target(im+(n-1)*n_DeltaM,:) = a*aratios(n,1)*[cos(F+DeltaM_vals(im));sin(F+DeltaM_vals(im));0];
            end
        end

        %RVrec = zeros(numel(t_vals),6,numel(alpha_vals));
        for l = 1:size(RV_out,2)

            for j = 1:numel(alpha_vals)
                %Integration
                alpha = alpha_vals(j);
                delta = delta_vals(j);
                for m = 1:n_coast
                    RV = RV_out(:,l);
                    for ti = 1:numel(t_vals)
                        t = t_vals(ti);
                        a_sail = asail(alpha,delta,RV);
                        a_sailxyz = projxyz(RV)*a_sail;
                        a_sail_coast = asail(alpha_coast(m),delta_coast(m),RV);
                        a_sail_coastxyz = projxyz(RV)*a_sail_coast;
                        if t_ahead*T < sail_max*T_out(l)
                            RV=rk4orbit(t,RV,tstep,a_sailxyz);
                        elseif t-tstart > t_ahead*T-sail_max*T_out(l)
                            RV=rk4orbit(t,RV,tstep,a_sailxyz);
                        else
                            RV=rk4orbit(t,RV,tstep,a_sail_coastxyz);
                        end
                    end
                    Rfinal1(:,2*j-1:2*j,l,m) = [RV(1:3)-[0;0;1e6],RV(1:3)+[0;0;1e6]];
                end
            end

            in = cell(size(RV_out,2),n_coast);
            F1 = cell(size(RV_out,2),n_coast);
            av1 = cell(size(RV_out,2),n_coast);

            for m = 1:n_coast
                try
                [F1{l,m},av1{l,m}] = convhull(Rfinal1(1,:,l,m),Rfinal1(2,:,l,m),Rfinal1(3,:,l,m));
                in{l,m} = CheckPointPolyhedron(Rfinal1(:,:,l,m)',F1{l,m},R_pla_target);
                catch
                    F1{l,m} = zeros(n_DeltaM*size(aratios,1),1);
                    av1{l,m} = zeros(n_DeltaM*size(aratios,1),1);
                    in{l,m} = zeros(n_DeltaM*size(aratios,1),1);
                end
                for im = 1:n_DeltaM*size(aratios,1)
                    if in{l,m}(im)
                        CATcell{iv}{mod(im-1,n_DeltaM)+1,floor((im-1)/n_DeltaM)+1} = [CATcell{iv}{mod(im-1,n_DeltaM)+1,floor((im-1)/n_DeltaM)+1}(:,:);[l,t_ahead]];
                    end
                    if ~caught(im) && in{l,m}(im)
                        caught(im) = 1;
                    end
                end
            end
        end

        if doPLOT
            clf;
            hold on;
            for l = 1:size(RV_out,2)
                for m = 1:n_coast
                    trisurf(F1{l,m}, Rfinal1(1,:,l,m), Rfinal1(2,:,l,m), Rfinal1(3,:,l,m),'FaceColor','cyan','FaceAlpha',0.2,'EdgeAlpha',0);
                    scatter3(Rfinal1(1,:,l,m), Rfinal1(2,:,l,m), Rfinal1(3,:,l,m), 10, 'filled', 'm');
                end
            end
            for i = 1:5
                plot3(R_path_pla{i}(:,1), R_path_pla{i}(:,2), R_path_pla{i}(:,3), 'm', 'LineWidth',0.5);
            end
            scatter3(0, 0, 0, 50, 'filled', 'MarkerFaceColor', 'y');

            scatter3(R_pla(:,1), R_pla(:,2), R_pla(:,3), 50, 'filled', 'm');
            scatter3(RV_flyby(1,1), RV_flyby(2,1), RV_flyby(3,1), 30, 'filled', 'b');
            scatter3(R_pla_target(:,1), R_pla_target(:,2), R_pla_target(:,3), 30, 'filled', 'b');
        end

        if doPLOT
            for im = 1:n_DeltaM*size(aratios,1)
                if caught(im) == 1
                    scatter3(R_pla_target(im,1), R_pla_target(im,2), R_pla_target(im,3), 30, 'filled', 'r');
                end
            end
        end

        if doPLOT
            hold off;

            axis equal;
            grid on;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title(sprintf('Stage 1: t\\_ahead = %.2f', t_ahead));
            xlim([-plotrange1,plotrange1]);
            ylim([-plotrange1,plotrange1]);
            zlim([-plotrange1,plotrange1]);
            drawnow;
        end
        fprintf('v in  = %1.2f, t = %1.3f Years\n',v_in/v_o,t_ahead);
    end
end

toc1 = toc;
fprintf('v_in_vals = %1.0f, N_t_vals = %1.0f, α-δ pairs = %1.0f, n_gamma = %1.0f, n_coast = %1.0f, RK steps = %1.0f\n',numel(v_in_vals),numel(t_ahead_vals),numel(alpha_vals),n_gamma,n_coast,RK_steps);
RKtot = numel(v_in_vals)*numel(alpha_vals)*numel(t_ahead_vals)*n_gamma*n_coast*RK_steps;
fprintf('Total RK steps = %1.0f\n',RKtot);
fprintf('Time spent = %1.2fs (thousand RK/sec = %1.5fs)',toc1,RKtot/10^3/toc1);

CATaxis = struct('vi_vals',v_in_vals/v_o,...
    'DeltaM_vals',DeltaM_vals,...
    'gamma_vals',gamma_vals,...
    'alpha_coast',alpha_coast(1:n_coast),...
    'delta_coast',delta_coast(1:n_coast),...
    'aratios',aratios...
    );

fname = strcat('CAT',string(datetime('now','Format',"yyyy-MM-dd-HH-mm-ss")),'.mat') ;
save(fname,'CATaxis','CATcell');

delete(gcp('nocreate'));   % stop all workers

doVIZ = 1;

if doVIZ


end

%% FUNCTIONS

function X = rk4orbit(t,X,h,a_sailxyz)
% Numerical Solution, Runge-Kutta 4th Order
k1 = forbit(t,X,a_sailxyz);
k2 = forbit(t+h/2,X+k1*h/2,a_sailxyz);
k3 = forbit(t+h/2,X+k2*h/2,a_sailxyz);
k4 = forbit(t+h,X+k3*h,a_sailxyz);

% Step forward in time
X = X+(h/6)*(k1+2*k2+2*k3+k4);

end

function Xdot = forbit(t,X,a_sailxyz)

mu = 139348062043.343;    % Gravitational constant (m^3/s^2)
r = X(1:3);                % Position (m)
v = X(4:6);                % Velocity (ms^2)

dr = v;
dv = (-mu/(norm(r))^3).*r + a_sailxyz; % Newton's law of gravity
Xdot = [dr; dv];

end

function a_sail = asail(alpha,delta,RV)
m = 500;
A = 15000;
C = 5.4026 * 10^-6;
r0 = 149597870.691; %[km]
r = norm(RV(1:3));
u_n = [cos(alpha),-sin(alpha),0;
    sin(alpha),cos(alpha),0;
    0,0,1]*[cos(-delta),0,sin(-delta);
    0,1,0;
    -sin(-delta),0,cos(-delta)]*[1,0,0]';
a_sail = 2*C*A/m*(r0/r)^2 * dot(u_n,[1,0,0])^2 * u_n./norm(u_n)/1000;
end

function proj = projxyz(RV)
rhat = RV(1:3)/norm(RV(1:3));
vhat = RV(4:6)/norm(RV(4:6));
nhat = cross(rhat,vhat)/norm(cross(rhat,vhat));
that = cross(nhat,rhat)/norm(cross(nhat,rhat));
proj = [rhat,that,nhat];
end

function I=CheckPointPolyhedron(V,F,vp)

for i=1:1:size(vp,1)

    omega=0;
    for j=1:1:size(F,1)
        omega=omega+Solid_Angle_Triangle(vp(i,:),V(F(j,1),:),V(F(j,2),:),V(F(j,3),:));
    end

    if abs(omega)>2*pi
        I(i)=true;
    else
        I(i)=false;
    end

end

end

function omega=Solid_Angle_Triangle(p,va,vb,vc)

r_pa=p-va; norm_r_pa=norm(r_pa);
r_pb=p-vb; norm_r_pb=norm(r_pb);
r_pc=p-vc; norm_r_pc=norm(r_pc);

omega=2*atan2(dot(r_pa,cross(r_pb,r_pc)),norm_r_pa*norm_r_pb*norm_r_pc+dot(r_pa,r_pb)*norm_r_pc+dot(r_pa,r_pc)*norm_r_pb+dot(r_pb,r_pc)*norm_r_pa);

end

function R = PlaPos(t,data,DeltaM)
a_pla = data(:,5);
ecc_pla = 0;
inc_pla     = 0/180*pi;
RAAN_pla       = 0/180*pi;
omega_pla  = 0/180*pi;
M0_pla = DeltaM;

mu = 139348062043.343; %[km^3/s^2]

n = sqrt(mu ./ a_pla.^3);

M = M0_pla(:) + n(:).*t;
E = M;
f_pla = 2*atan(sqrt((1+ecc_pla)./(1-ecc_pla)).*tan(E/2));

r_pla = a_pla.*(1-ecc_pla.^2)./(1+ecc_pla.*cos(f_pla));

R = [r_pla.*(cos(f_pla+omega_pla).*cos(RAAN_pla)-sin(f_pla+omega_pla).*cos(inc_pla).*sin(RAAN_pla)),...
    r_pla.*(cos(f_pla+omega_pla).*sin(RAAN_pla)+sin(f_pla+omega_pla).*cos(inc_pla).*cos(RAAN_pla)),...
    r_pla.*(sin(f_pla+omega_pla).*sin(inc_pla))];


end

function [arec,erec,f0rec] = lambert(r1,r2,eps,mu,TOF);

tol = 0.01;
sgn = sign(r2 - r1);
inv_mu = 1/sqrt(mu);
n_loops_max = 10;

f0min = atan((cos(eps)-r1/r2)/sin(eps));
A = (r1/r2) - cos(eps);
B = sin(eps);
f0max1 = atan2(B, A) + acos((1 - (r1/r2))/sqrt(A^2 + B^2));
f0max2 = atan2(B, A) - acos((1 - (r1/r2))/sqrt(A^2 + B^2));
f0maxR = atan((cos(eps)-1)/sin(eps));

if eps < pi
    if r2 > r1
        f0range = [f0min,f0max1];
        f0rangeMini = [f0max1,f0max2];
    else
        f0range = [f0max1-2*pi,f0min];
        f0rangeMini = [f0max1-2*pi,f0max2];
    end
else
    if r2 > r1
        f0range = [f0maxR-pi,f0max1];
        f0rangeMini = [f0max1,f0max2];
    else
        f0range = [f0max1-2*pi,f0maxR-pi];
        f0rangeMini = [f0max1-2*pi,f0max2];
    end
end

n_calc = 0;

while true
    n_calc = n_calc+1;
    f0 = (f0range(1)+f0range(2))/2;
    ecc = ((1-r1/r2)/(cos(f0)*r1/r2-cos(f0+eps)));
    a = r1*(1+ecc*cos(f0))/(1-ecc^2);
    if ecc < 1
        k = sqrt((1 - ecc) / (1 + ecc));
        E0 = 2 * atan(k * tan(f0/2));
        E2 = 2 * atan(k * tan((f0 + eps)/2));
        M0 = (E0-ecc*sin(E0));
        M2 = (E2-ecc*sin(E2));
    else
        k = sqrt((ecc - 1) / (ecc +1));
        H0 = 2 * atanh(k * tan(f0/2));
        H2 = 2 * atanh(k * tan((f0 + eps)/2));
        M0 = (ecc*sinh(H0)-H0);
        M2 = (ecc*sinh(H2)-H2);
    end
    M = M2-M0 + 2*pi*(M2<M0);
    tof = M*abs(a)*sqrt(abs(a)/mu);
    if abs(tof - TOF) < tol
        break;
    end
    if sgn*(tof - TOF) > 0
        f0range(2) = f0;
    else
        f0range(1) = f0;
    end

end

arec = a;
erec = ecc;
f0rec = f0;

endloop = 0;
skip = 0;
for loop = 1:n_loops_max
    if skip
        break;
    end

    f0rangeM = f0rangeMini;
    n = 0;
    while true
        n = n+1;
        n_calc = n_calc+1;
        f0(1) = (f0rangeM(1)+f0rangeM(2))/2 - 0.0001;
        f0(2) = (f0rangeM(1)+f0rangeM(2))/2 + 0.0001;
        ecc = ((1-r1/r2)./(cos(f0)*r1/r2-cos(f0+eps)));
        a = r1*(1+ecc.*cos(f0))./(1-ecc.^2);
        k = sqrt((1 - ecc) ./ (1 + ecc));
        E0 = 2 * atan(k .* tan(f0/2));
        E2 = 2 * atan(k .* tan((f0 + eps)/2));
        M0 = (E0-ecc.*sin(E0));
        M2 = (E2-ecc.*sin(E2));
        M = M2-M0 + 2*pi*(M2(1)<M0(1));
        tof = (2*pi*loop+M)*a(2)*sqrt(a(2))*inv_mu;
        if tof(1) < TOF
            break;
        else
            if tof(2) > tof(1)
                f0rangeM(1) = f0(1);
            else
                f0rangeM(2) = f0(1);
            end
        end
        if n > 12
            endloop = 1;
            break;
        end
    end

    if endloop
        break;
    end

    f0min = f0(1);
    tofmin = tof(1);

    f0rangeMleft = [f0rangeM(1),f0min];
    f0rangeMright = [f0min,f0rangeM(2)];
    while true
        n_calc = n_calc+2;
        f0(1) = (f0rangeMleft(1)+f0rangeMleft(2))/2;
        f0(2) = (f0rangeMright(1)+f0rangeMright(2))/2;
        ecc = ((1-r1/r2)./(cos(f0)*r1/r2-cos(f0+eps)));
        a = r1*(1+ecc.*cos(f0))./(1-ecc.^2);
        k = sqrt((1 - ecc) ./ (1 + ecc));
        E0 = 2 * atan(k .* tan(f0/2));
        E2 = 2 * atan(k .* tan((f0 + eps)/2));
        M0 = (E0-ecc.*sin(E0));
        M2 = (E2-ecc.*sin(E2));
        M = M2-M0 + 2*pi*(M2(1)<M0(1));
        tof = (2*pi*[loop,loop-1]+M).*(a).*sqrt(a)*inv_mu;
        if tof(1) > TOF
            f0rangeMleft(1) = f0(1);
        else
            f0rangeMleft(2) = f0(1);
        end
        if tof(2) > TOF
            f0rangeMright(2) = f0(2);
        else
            f0rangeMright(1) = f0(2);
        end
        if abs(tof(1) - TOF) < tol && abs(tof(2) - TOF) < tol
            break;
        end
        if n_calc > 1000
            %fprintf(';-; %1.f %1.f\n',r2,loop);
            skip = 1;
            break;
        end
    end
    arec = [arec,a];
    erec = [erec,ecc];
    f0rec = [f0rec,f0];
end
end

function COE = RV2COE(RV, mu)
%RV2COE Converts state vectors (r,v) to classical orbital elements
%   RV: 6xn matrix [r; v]
%   mu: gravitational parameter
%   COE: [a; e; i; RAAN; omega; M0] (6xn)

n = size(RV, 2);
COE = zeros(6, n);

for k = 1:n
    r = RV(1:3, k);
    v = RV(4:6, k);

    R = norm(r);
    V = norm(v);

    h = cross(r, v);
    h_mag = norm(h);

    k_hat = [0 0 1]';
    n_vec = cross(k_hat, h);
    n_mag = norm(n_vec);

    e_vec = (1/mu) * ((V^2 - mu/R) * r - dot(r, v) * v);
    e = norm(e_vec);

    E = V^2/2 - mu/R;
    if abs(E) > 1e-12
        a = -mu / (2*E);
    else
        a = inf; % parabolic case
    end

    i = acos(h(3)/h_mag);

    if n_mag ~= 0
        RAAN = atan2(n_vec(2), n_vec(1));
    else
        RAAN = 0;
    end

    if n_mag ~= 0 && e > 1e-10
        omega = atan2(dot(cross(n_vec, e_vec), h)/h_mag, dot(n_vec, e_vec));
    else
        omega = 0;
    end

    if e > 1e-10
        f = atan2(dot(cross(e_vec, r), h)/(h_mag*e), dot(e_vec, r)/(e*R));
    else
        f = atan2(r(2), r(1));
    end

    Ecc = 2*atan( sqrt((1-e)/(1+e)) * tan(f/2) );
    M0 = Ecc - e*sin(Ecc);

    COE(:,k) = [a; e; i; RAAN; omega+pi; M0];
end
end

function R_path = COE2Path(COE, f_vals, mu)
%COE2Path Converts COE(6xn) to orbital position paths
%   COE: [a; e; i; RAAN; omega; M0] (6xn)
%   f_vals: vector of true anomaly values (e.g. 0:0.01:2*pi)
%   mu: gravitational parameter
%   R_path: [length(f_vals) x 3 x n] position vectors

n = size(COE, 2);
R_path = zeros(length(f_vals), 3, n);

for k = 1:n
    a = COE(1,k);
    e = COE(2,k);
    inc = COE(3,k);
    RAAN = COE(4,k);
    omega = COE(5,k)-pi;
    r_path = a.*(1-e.^2)./(1+e.*cos(f_vals));
    R_path(:,:,k) = [r_path.*(cos(f_vals+omega)*cos(RAAN)-sin(f_vals+omega)*cos(inc)*sin(RAAN)),...
        r_path.*(cos(f_vals+omega)*sin(RAAN)+sin(f_vals+omega)*cos(inc)*cos(RAAN)),...
        r_path.*(sin(f_vals+omega)*sin(inc))];
end
end



