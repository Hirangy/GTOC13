clear; clc;

% === Load your result file ===
load('Results02-47336400d0/ResultP02-47336400d0-P02-66750623d728.mat'); 
% expect variable: result_table  (N x 12)

if ~exist('result_table','var')
    error('result_table not found in file.');
end

% Columns: 1:body_id, 2:flag, 3:t, 4-6:r, 7-9:v, 10-12:u
t = result_table(:,3);        % Nx1
r = result_table(:,4:6);      % Nx3
v = result_table(:,7:9);      % Nx3
u = result_table(:,10:12);    % Nx3

mu = 139348062043.343; % m^3/s^2 (keeps your constant)
% sail parameters (used only for magnitude scaling if needed)
A = 15000; m = 500; C = 5.4026e-6;
r0 = 149597870.691; % km

fprintf('Checking RK4 error tolerance (1e-4) — fixed indexing...\n\n');
tol = 1e-4;
pass_count = 0;
seg_count = 0;

for i = 1:(size(result_table,1)-1)
    t0 = t(i);
    t1 = t(i+1);

    % skip identical-time rows (control discontinuities)
    if t1 == t0
        continue;
    end

    seg_count = seg_count + 1;

    X0 = [r(i,:)'; v(i,:)'];     % 6x1
    Xf = [r(i+1,:)'; v(i+1,:)']; % 6x1

    dt = t1 - t0;
    substeps = 100;
    h = dt / substeps;

    u0 = u(i,:)';
    u1 = u(i+1,:)';

    % protect against zero vectors
    if norm(u0) == 0 || norm(u1) == 0
        warning('Zero u vector at segment %d — skipping', i);
        continue;
    end

    u0 = u0 / norm(u0);
    u1 = u1 / norm(u1);

    X = X0;

    for k = 1:substeps
        frac = (k-1) / (substeps-1); % goes 0 .. 1 across substeps
        u_interp = (1 - frac) * u0 + frac * u1;
        u_interp = u_interp / norm(u_interp);

        % Use u_interp as the sail normal direction. Need only magnitude if
        % your asail formula expects angles; here we compute acceleration
        % magnitude using your asail formula reduced to u-direction:
        rmag = norm(X(1:3));
        % compute scalar factor from asail (assuming dot(u,body_x)=cosTheta)
        cosTheta = dot(u_interp, [1;0;0]); % body x assumed [1;0;0] in asail deriv
        a_mag = 2*C*A/m * (r0 / rmag)^2 * (cosTheta)^2; % in km/s^2? asail had /1000
        a_sailxyz = a_mag * u_interp / 1000; % match your asail (/1000 conversion)

        % integrate one RK4 step using your rk4orbit signature
        X = rk4orbit(t0 + (k-1)*h, X, h, a_sailxyz);
    end

    X_RK4 = X;

    denom = norm(Xf - X0);
    if denom == 0
        err_rel = norm(X_RK4 - Xf); % avoid div0, absolute error
    else
        err_rel = norm(X_RK4 - Xf) / denom;
    end

    pass = err_rel < tol;
    if pass
        pass_count = pass_count + 1;
    end

    fprintf('Segment %3d (t=%.3f -> %.3f s): rel err = %.3e %s\n', ...
        seg_count, t0, t1, err_rel, tern(pass,'✅ PASS','❌ FAIL'));
end

fprintf('\nSummary: %d/%d segments passed (%.1f%%)\n', ...
    pass_count, seg_count, 100*pass_count/max(1,seg_count));


% --- small helper inline function for ternary ----------
function s = tern(cond, a, b)
    if cond, s = a; else s = b; end
end


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