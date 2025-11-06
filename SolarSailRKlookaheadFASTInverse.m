clear;
clc;

data = readmatrix('/Users/hirangy/Downloads/gtoc13_asteroids.csv', 'NumHeaderLines', 1);
data = data(1:257,:);

%INITIALIZATION
doPLOT = 1;
mu = 139348062043.343; %[km^3/s^2]
r0 = 149597870.691; %[km]
a = 3.8*r0;
e = 0;
inc = 0;
RAAN = 2;
omega = 0;
M0 = 2;

%Planned duration of sailing (fraction of SC Period)
sail_dur = 0.8;

%Original SC path

f_vals = (0:0.001:2*pi)';
r_path = a.*(1-e.^2)./(1+e.*cos(f_vals));
R_path = [r_path.*(cos(f_vals+omega)*cos(RAAN)-sin(f_vals+omega)*cos(inc)*sin(RAAN)),...
    r_path.*(cos(f_vals+omega)*sin(RAAN)+sin(f_vals+omega)*cos(inc)*cos(RAAN)),...
    r_path.*(sin(f_vals+omega)*sin(inc))];

T = 2*pi*sqrt(a^3/mu);
n = sqrt(mu/a^3);
p = a*(1-e^2);


%Initial net with 9 points
delta_vals = [0,50/180*pi,25/180*pi,-50/180*pi,-25/180*pi];
delta_vals = [delta_vals,delta_vals,delta_vals,delta_vals,0,0,pi/2];
alpha_vals = [zeros(1,5)+0.5,zeros(1,5)-0.25,zeros(1,5)+0.2,zeros(1,5)-0.25,-1,1,0];

M = M0;
E = M;
for iter = 1:30
    E = E - (E - e.*sin(E) - M)./(1 - e.*cos(E));
end
f = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
r = a.*(1-e.^2)./(1+e.*cos(f));
v = sqrt(2*mu./r-mu/a);

%Original RV for keepsake
RVori(1:3,1) = [r.*(cos(f+omega)*cos(RAAN)-sin(f+omega)*cos(inc)*sin(RAAN)),...
    r.*(cos(f+omega)*sin(RAAN)+sin(f+omega)*cos(inc)*cos(RAAN)),...
    r.*(sin(f+omega)*sin(inc))];
gamma = atan(e.*sin(f)./(1+e.*cos(f)));
RVori(4:6,1) = [v.*(-sin(f+omega-gamma)*cos(RAAN)-cos(f+omega-gamma)*cos(inc)*sin(RAAN)),...
    v.*(-sin(f+omega-gamma)*sin(RAAN)+cos(f+omega-gamma)*cos(inc)*cos(RAAN)),...
    v.*(cos(f+omega-gamma)*sin(inc))];
Hori = cross(RVori(1:3),RVori(4:6))/norm(cross(RVori(1:3),RVori(4:6)));

plotrange1 = 5*r0;
tic;
%FIRST STAGE: sweep large arc to catch any asteroid
if doPLOT
    figure('Position', [100 100 1100 700]);
end
t_ahead_vals = sail_dur:0.01:2;
for k = 1:numel(t_ahead_vals)
    t_ahead = t_ahead_vals(k);
    tstep = t_ahead*T/1000;
    t_vals = 0:tstep:t_ahead*T;
    R_ast = AstPos(t_ahead*T,data);
    Rfinal1 = zeros(3,numel(alpha_vals));
    %RVrec = zeros(numel(t_vals),6,numel(alpha_vals));
    for j = 1:numel(alpha_vals)
        RV = RVori;
        %Integration
        alpha = alpha_vals(j);
        delta = delta_vals(j);
        for ti = 1:numel(t_vals)
            t = t_vals(ti);
            a_sail = asail(alpha,delta,RV);
            a_sailxyz = projxyz(RV)*a_sail;
            if t > (t_ahead-sail_dur)*T
                RV=rk4orbit(t,RV,tstep,a_sailxyz);
            else
                RV=rk4orbit(t,RV,tstep,[0,0,0]');
            end
            %RVrec(ti,:,j) = RV';
        end
        Rfinal1(:,j) = RV(1:3);
    end
    [F1,av1] = convhull(Rfinal1(1,:),Rfinal1(2,:),Rfinal1(3,:));
    in = CheckPointPolyhedron(Rfinal1',F1,R_ast);

    if doPLOT
        clf;
        hold on;
        plot3(R_path(:,1), R_path(:,2), R_path(:,3), 'b');
        scatter3(0, 0, 0, 50, 'filled', 'MarkerFaceColor', 'y');
        trisurf(F1, Rfinal1(1,:), Rfinal1(2,:), Rfinal1(3,:), ...
            'FaceColor','cyan','FaceAlpha',0.2,'EdgeAlpha',0);
        scatter3(R_ast(:,1), R_ast(:,2), R_ast(:,3), 20, 'filled', 'g');
        scatter3(RVori(1),RVori(2),RVori(3), 50, 'filled', 'MarkerFaceColor', 'b');
    end
    if any(in)
        if doPLOT
            scatter3(R_ast(in,1), R_ast(in,2), R_ast(in,3), 30, 'filled', 'r');
        end
        i_ast = find(in,1);
        R_ast_in = R_ast(in,:);
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

    if any(in) 
        text(R_ast_in(1),R_ast_in(2),R_ast_in(3),'Found one!', 'Fontsize',20,'HorizontalAlignment','center');
        break;
    end
end

if ~any(in)
hold on;
text(0.2,0.5,'DIDNT FIND NO ASTEROID (╥﹏╥) huh', 'Fontsize',20,'HorizontalAlignment','center');
drawnow;
pause(5);
end
t_ahead_break = t_ahead;
toc1 = toc;
tic;
%Alpha and delta net now refined to 9x13 values

delta_vals = [0,10,25,40,50,75,-10,-25,-40,-50,-75]/180*pi;
alpha_vals = [pi/2,0,0.125,0.25,0.5,0.75,1,1.25,1.4,-0.125,-0.25,-0.5,-0.75,-1,-1.25,1.4];
[Alpha, Delta] = meshgrid(alpha_vals, delta_vals);
alpha_vals = Alpha(:);
delta_vals = Delta(:);

%SECOND STAGE: see where the asteroid touches the bubble

if doPLOT
    figure('Position', [100 100 1100 700]);
end
plotrange2 = 0.5*r0;
t_ahead_vals = [t_ahead_break:-0.003:t_ahead_break-0.015,t_ahead_break-0.020:-0.005:t_ahead_break-0.050];
nk = numel(t_ahead_vals);
for k = 1:numel(t_ahead_vals)
    t_ahead = t_ahead_vals(k);
    t_vals = linspace(0,t_ahead*T,500);
    tstep = t_ahead*T/500;
    R_ast = AstPos(t_ahead*T,data);
    Rfinal1 = zeros(3,numel(alpha_vals));
    %RVrec = zeros(numel(t_vals),6,numel(alpha_vals));
    for j = 1:numel(alpha_vals)
        RV = RVori;
        %Integration
        alpha = alpha_vals(j);
        delta = delta_vals(j);
        for ti = 1:numel(t_vals)
            t = t_vals(ti);
            a_sail = asail(alpha,delta,RV);
            a_sailxyz = projxyz(RV)*a_sail;
            if t > (t_ahead-sail_dur)*T
                RV=rk4orbit(t,RV,tstep,a_sailxyz);
            else
                RV=rk4orbit(t,RV,tstep,[0,0,0]');
            end
            %RVrec(ti,:,j) = RV';
        end
        Rfinal1(:,j) = RV(1:3);
    end
    [F1,av1] = convhull(Rfinal1(1,:),Rfinal1(2,:),Rfinal1(3,:));
    in = CheckPointPolyhedron(Rfinal1',F1,R_ast(i_ast,:));
    if doPLOT
        clf;  % clear current figure
        hold on;
        plot3(R_path(:,1), R_path(:,2), R_path(:,3), 'b');
        scatter3(0, 0, 0, 50, 'filled', 'MarkerFaceColor', 'y');
        trisurf(F1, Rfinal1(1,:), Rfinal1(2,:), Rfinal1(3,:), ...
            'FaceColor','cyan','FaceAlpha',0.2,'EdgeAlpha',0);
        scatter3(R_ast(:,1), R_ast(:,2), R_ast(:,3), 20, 'filled', 'g');
    end
    R_ast_in = R_ast(i_ast,:);

    if ~in
        if doPLOT
            scatter3(R_ast(i_ast,1), R_ast(i_ast,2), R_ast(i_ast,3), 30, 'filled', 'r');
        end
    end
    for ia = 1:numel(alpha_vals)
        dist(ia) =norm(R_ast_in'-Rfinal1(:,ia));
    end
    [dist_min,imin] = min(dist);

    if doPLOT
        scatter3(Rfinal1(1,:), Rfinal1(2,:), Rfinal1(3,:), 10,'b','filled');
        scatter3(R_ast(i_ast,1), R_ast(i_ast,2), R_ast(i_ast,3), 30, 'filled', 'r');
    end
    %SECOND LOOP SUBSAMPLING: take α-δ pair and subdivide
    if ~in
        delta_vals = linspace(delta_vals(imin)-pi/10,delta_vals(imin)+pi/12,10);
        alpha_vals = linspace(alpha_vals(imin)-pi/10,alpha_vals(imin)+pi/12,10);
        [Alpha, Delta] = meshgrid(alpha_vals, delta_vals);
        alpha_vals = Alpha(:);
        delta_vals = Delta(:);
        Rfinal1 = zeros(3,numel(alpha_vals));
        for j = 1:numel(alpha_vals)
            RV = RVori;
            %Integration
            alpha = alpha_vals(j);
            delta = delta_vals(j);
            for ti = 1:numel(t_vals)
                t = t_vals(ti);
                a_sail = asail(alpha,delta,RV);
                a_sailxyz = projxyz(RV)*a_sail;
                if t > (t_ahead-sail_dur)*T
                    RV=rk4orbit(t,RV,tstep,a_sailxyz);
                else
                    RV=rk4orbit(t,RV,tstep,[0,0,0]');
                end
                %RVrec(ti,:,j) = RV';
            end
            Rfinal1(:,j) = RV(1:3);
        end
        if doPLOT
            scatter3(Rfinal1(1,:), Rfinal1(2,:), Rfinal1(3,:), 10 ,'b');
        end
        for ia = 1:numel(Rfinal1)/3
            dist(ia) =norm(R_ast_in'-Rfinal1(:,ia));
        end
        [dist_min,imin] = min(dist);
        if doPLOT
            plot3([R_ast(i_ast,1),Rfinal1(1,imin)],[R_ast(i_ast,2),Rfinal1(2,imin)],[R_ast(i_ast,3),Rfinal1(3,imin)]);
        end
    end

    if doPLOT
        hold off;
        axis equal;
        grid on;
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(sprintf('Stage 2: t\\_ahead = %.3f', t_ahead));
        xlim(R_ast_in(1)+[-plotrange2,plotrange2]);
        ylim(R_ast_in(2)+[-plotrange2,plotrange2]);
        zlim(R_ast_in(3)+[-plotrange2,plotrange2]);
        drawnow;
    end
    if ~in
        break;
    end
end
toc2 = toc;
tic;

%{
Timewise_step = [1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12];
Timewise_range = [30,15,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12];
Timewise_size = 4*[1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9];
Timewise_res = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5];
Spacewise_size = 14*[1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9];
Spacewise_res = [8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8];
%}

%THIRD STAGE: subdividing the surface until locally linear
base = 6;
Timewise_step = 5*1e-4*base.^(0:-1:-30);
Timewise_range = [10,5,10,10,10,10,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8];
Timewise_size = 20*base.^(0:1:30);
Timewise_res = 3+zeros(1,30);
Spacewise_size = 3*Timewise_size;
Spacewise_res = 5+zeros(1,30);

for loop = 1:1
    delta_vals = [linspace(delta_vals(imin)-pi/Timewise_size(loop),delta_vals(imin)+pi/Timewise_size(loop),Timewise_res(loop)),50/180*pi,25/180*pi,-50/180*pi,-25/180*pi];
    alpha_vals = [linspace(alpha_vals(imin)-pi/Timewise_size(loop),alpha_vals(imin)+pi/Timewise_size(loop),Timewise_res(loop)),0.5,-0.5,0.25,-0.25];
    [Alpha, Delta] = meshgrid(alpha_vals, delta_vals);
    alpha_vals = Alpha(:);
    delta_vals = Delta(:);
    alpha_vals = [alpha_vals;0];
    delta_vals = [delta_vals;pi/2];

    %THIRD (3-max) LOOP: SAME AS SECOND LOOP JUST SMALLER AND MORE PRECISE
    t_ahead_break = t_ahead;
    if doPLOT
    figure('Position', [100 100 1100 700]);
    end
    plotrange2 = 2/Timewise_size(loop)*r0;
    t_ahead_vals = t_ahead_break-Timewise_step(loop)*Timewise_range(loop):Timewise_step(loop):t_ahead_break+Timewise_step(loop)*Timewise_range(loop)+Timewise_step(loop)*Timewise_range(loop)/2;
    nk = numel(t_ahead_vals);
    clear Rfinal1;
    for k = 1:numel(t_ahead_vals)
        t_ahead = t_ahead_vals(k);
        t_vals = linspace(0,t_ahead*T,1000);
        tstep = t_ahead*T/1000;
        R_ast = AstPos(t_ahead*T,data);
        Rfinal1 = zeros(3,numel(alpha_vals));
        %RVrec = zeros(numel(t_vals),6,numel(alpha_vals));
        for j = 1:numel(alpha_vals)
            RV = RVori;
            %Integration
            alpha = alpha_vals(j);
            delta = delta_vals(j);
            for ti = 1:numel(t_vals)
                t = t_vals(ti);
                a_sail = asail(alpha,delta,RV);
                a_sailxyz = projxyz(RV)*a_sail;
                if t > (t_ahead-sail_dur)*T
                    RV=rk4orbit(t,RV,tstep,a_sailxyz);
                else
                    RV=rk4orbit(t,RV,tstep,[0,0,0]');
                end
                %RVrec(ti,:,j) = RV';
            end
            Rfinal1(:,j) = RV(1:3);
        end
        [F1,av1] = convhull(Rfinal1(1,:),Rfinal1(2,:),Rfinal1(3,:));
        in = CheckPointPolyhedron(Rfinal1',F1,R_ast(i_ast,:));
        clf;  % clear current figure
        hold on;
        plot3(R_path(:,1), R_path(:,2), R_path(:,3), 'b');
        scatter3(0, 0, 0, 50, 'filled', 'MarkerFaceColor', 'y');
        trisurf(F1, Rfinal1(1,:), Rfinal1(2,:), Rfinal1(3,:), ...
            'FaceColor','cyan','FaceAlpha',0.2,'EdgeAlpha',0);
        scatter3(R_ast(:,1), R_ast(:,2), R_ast(:,3), 20, 'filled', 'g');
        scatter3(Rfinal1(1,:), Rfinal1(2,:), Rfinal1(3,:), 10 ,'b');
        R_ast_in = R_ast(i_ast,:);
        scatter3(R_ast(i_ast,1), R_ast(i_ast,2), R_ast(i_ast,3), 30, 'filled', 'r');

        clear dist;
        for ia = 1:numel(alpha_vals)
            dist(ia) =norm(R_ast_in'-Rfinal1(:,ia));
        end
        [dist_min,imin] = min(dist);

        %THIRD LOOP SUBSAMPLING: take α-δ pair and subdivide
        if in
            if doPLOT
                scatter3(R_ast(i_ast,1), R_ast(i_ast,2), R_ast(i_ast,3), 30, 'filled', 'r');
            end
            delta_vals = linspace(delta_vals(imin)-pi/Spacewise_size(loop),delta_vals(imin)+pi/Spacewise_size(loop),Spacewise_res(loop));
            alpha_vals = linspace(alpha_vals(imin)-pi/Spacewise_size(loop),alpha_vals(imin)+pi/Spacewise_size(loop),Spacewise_res(loop));
            [Alpha, Delta] = meshgrid(alpha_vals, delta_vals);
            alpha_vals = Alpha(:);
            delta_vals = Delta(:);
            Rfinal1 = zeros(3,numel(alpha_vals));
            for j = 1:numel(alpha_vals)
                RV = RVori;
                %Integration
                alpha = alpha_vals(j);
                delta = delta_vals(j);
                for ti = 1:numel(t_vals)
                    t = t_vals(ti);
                    a_sail = asail(alpha,delta,RV);
                    a_sailxyz = projxyz(RV)*a_sail;
                    if t > (t_ahead-sail_dur)*T
                        RV=rk4orbit(t,RV,tstep,a_sailxyz);
                    else
                        RV=rk4orbit(t,RV,tstep,[0,0,0]');
                    end
                    %RVrec(ti,:,j) = RV';
                end
                Rfinal1(:,j) = RV(1:3);
            end
            if doPLOT
                scatter3(Rfinal1(1,:), Rfinal1(2,:), Rfinal1(3,:), 10 ,'b','filled');
            end
            clear dist;
            for ia = 1:numel(Rfinal1)/3
                dist(ia) = norm(R_ast_in'-Rfinal1(:,ia));
            end
            [dist_min,imin] = min(dist);
            if doPLOT
                plot3([R_ast(i_ast,1),Rfinal1(1,imin)],[R_ast(i_ast,2),Rfinal1(2,imin)],[R_ast(i_ast,3),Rfinal1(3,imin)]);
            end
        end
        if doPLOT

            hold off;
            axis equal;
            grid on;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title(sprintf('Stage 3: t\\_ahead = %.12f', t_ahead));
            xlim(R_ast_in(1)+[-plotrange2,plotrange2]);
            ylim(R_ast_in(2)+[-plotrange2,plotrange2]);
            zlim(R_ast_in(3)+[-plotrange2,plotrange2]);
            drawnow;
        end
        if in
            fprintf('Closest Trajectory = %1.3f km\n', dist_min);
            break;
        end
    end
end
toc3 = toc;
tic;
t_res_la = 10000;

%FORTH STAGE: using linear interpolation to reach convergence
if doPLOT
    figure('Position', [100 100 1100 700]);
    text(0.2,0.5,'LINEAR ALGEBRA IN PROGRESS (∪｡∪)｡｡｡zzZ	', 'Fontsize',20);
    title('Stage 4');
    hold off;
    drawnow;
end
Rf_min = Rfinal1(:,imin);
Delta_R = R_ast_in'-Rf_min;
alpha_min = alpha_vals(imin);
delta_min = delta_vals(imin);

for laloop = 1:15

alpha_vals=[alpha_min,alpha_min+1e-9,alpha_min,alpha_min];
delta_vals=[delta_min,delta_min,delta_min+1e-9,delta_min];
t_ahead_vals = [t_ahead,t_ahead,t_ahead,t_ahead+1e-9];
Rfinal1 = zeros(3,4);

for j = 1:numel(alpha_vals)
    t_vals = linspace(0,t_ahead_vals(j)*T,t_res_la);
    RV = RVori;
    %Integration
    alpha = alpha_vals(j);
    delta = delta_vals(j);
    tstep = t_ahead_vals(j)*T/t_res_la;
    for ti = 1:numel(t_vals)
        t = t_vals(ti);
        a_sail = asail(alpha,delta,RV);
        a_sailxyz = projxyz(RV)*a_sail;
        if t > (t_ahead-sail_dur)*T
            RV=rk4orbit(t,RV,tstep,a_sailxyz);
        else
            RV=rk4orbit(t,RV,tstep,[0,0,0]');
        end
        %RVrec(ti,:,j) = RV';
    end
    Rfinal1(1:3,j) = RV(1:3);
end

Ra = (Rfinal1(:,2)-Rfinal1(:,1))/1e-9;
Rd = (Rfinal1(:,3)-Rfinal1(:,1))/1e-9;
Rt = (Rfinal1(:,4)-Rfinal1(:,1))/1e-9;

J = [Ra Rd Rt];
dP = J \ Delta_R;
alpha_new = alpha_min + dP(1);
delta_new = delta_min + dP(2);
t_ahead_new = t_ahead + dP(3);

alpha_vals = alpha_new;
delta_vals = delta_new;
t_ahead = t_ahead_new;

Rfinal1 = zeros(1,3);
t_vals = linspace(0,t_ahead*T,t_res_la);
for j = 1:numel(alpha_vals)
    RV = RVori;
    %Integration
    alpha = alpha_vals(j);
    delta = delta_vals(j);
    tstep = t_ahead*T/t_res_la;
    for ti = 1:numel(t_vals)
        t = t_vals(ti);
        a_sail = asail(alpha,delta,RV);
        a_sailxyz = projxyz(RV)*a_sail;
        if t > (t_ahead-sail_dur)*T
            RV=rk4orbit(t,RV,tstep,a_sailxyz);
        else
            RV=rk4orbit(t,RV,tstep,[0,0,0]');
        end
        %RVrec(ti,:,j) = RV';
    end
    Rfinal1(1:3,j) = RV(1:3);
    dist = norm(R_ast_in'-Rfinal1(1:3,1));
    fprintf('Linear Interpolation %1.f = %1.3f km\n', laloop ,dist);

end
Delta_R = R_ast_in'-Rfinal1;
alpha_min = alpha;
delta_min = delta;
    if dist < 1e-5
        break;
    end
end

if doPLOT
    figure('Position', [100 100 1100 700]);
    text(0.2,0.5,'LINEAR ALGEBRA IN PROGRESS...', 'Fontsize',20);
    hold on;
    if dist > 0.1
    text(0.2,0.3,'did not converge (ﾉಥ益ಥ)ﾉ', 'Fontsize',20);
    else
    text(0.2,0.3,'DONE (o\^▽\^o)', 'Fontsize',20);
    end
    title('Stage 4');
end

if dist > 0.1
    fprintf('Final Close Approach = %1.3f km\n', dist);
elseif dist > 0.0001
    fprintf('Final Close Approach = %1.3f m\n', dist*1000);
elseif dist > 0.0000001
    fprintf('Final Close Approach = %1.3f mm\n', dist*1000000);
else
    fprintf('Final Close Approach = %1.3f μm\n', dist*1000000000);
end
toc4 = toc;
fprintf('Calculation time = %1.2f + %1.2f + %1.2f + %1.2f s = %1.2f s\n', toc1,toc2,toc3,toc4,toc1+toc2+toc3+toc4);
fprintf('α = %1.7f deg, δ = %1.7f deg, T_act = %1.7f d, T_sail = %1.2fT\n', alpha_new*180/pi, delta_new*180/pi, (t_ahead-sail_dur)*T/86400, sail_dur);


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
vhat = RV(4:6)/norm(RV(1:3));
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

function R = AstPos(t,data)

ast_id     = data(:,1);
a_ast = data(:,2);
ecc_ast  = data(:,3);
inc_ast      = data(:,4)/180*pi;
RAAN_ast       = data(:,5)/180*pi;
omega_ast  = data(:,6)/180*pi;
M0_ast = data(:,7)/180*pi;
weight_ast        = data(:,8);

mu = 139348062043.343; %[km^3/s^2]
AU = 149597870.691; %[km]

n = sqrt(mu ./ a_ast.^3);

M = M0_ast(:) + n(:).*t;
E = M;
for iter = 1:20
    E = E - (E - ecc_ast.*sin(E) - M)./(1 - ecc_ast.*cos(E));
end
f_ast = 2*atan(sqrt((1+ecc_ast)./(1-ecc_ast)).*tan(E/2));

r_ast = a_ast.*(1-ecc_ast.^2)./(1+ecc_ast.*cos(f_ast));

R = [r_ast.*(cos(f_ast+omega_ast).*cos(RAAN_ast)-sin(f_ast+omega_ast).*cos(inc_ast).*sin(RAAN_ast)),...
    r_ast.*(cos(f_ast+omega_ast).*sin(RAAN_ast)+sin(f_ast+omega_ast).*cos(inc_ast).*cos(RAAN_ast)),...
    r_ast.*(sin(f_ast+omega_ast).*sin(inc_ast))];


end
