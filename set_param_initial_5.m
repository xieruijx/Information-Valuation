%% General parameter

error = 1e-8;

Num_t = 24;

Num_m = 30;

m_C = 0.02;
m_sigma = 1; %
m_demand_mean = [0.44; 0.46; 0.51; 0.58; 0.60; 0.64;
    0.66; 0.70; 0.76; 0.80; 0.81; 0.84;
    0.80; 0.76; 0.72; 0.61; 0.60; 0.63;
    0.69; 0.80; 0.77; 0.64; 0.56; 0.51] / 0.84;


delta = 0.95; %
xi = delta; %

Gamma_T = sqrt(Num_t*(1-delta)*(1+Num_t-Num_t*xi)/(1-xi));

%% Parameter for generator

Num_i = 3;

alpha = [30;40;20]*ones(1,Num_t); % dollar/1 MW
gamma = alpha;

p_max = [5.2;2;6]; % MW
p_min = 0.3 * p_max;
r_max = 0.4 * p_max; % reserve
R_max = 0.4 * p_max; % ramp

%% Parameter for prosumer

Num_j = 5;
% m_demand_mean = normrnd(ones(Num_j,1)*m_demand_mean',0.01*ones(Num_j,Num_t));
% d_mean = ([0; -1; 3; 3; 4] * ones(1,Num_t)) .* m_demand_mean;
% save('d_mean.mat','d_mean');
load('d_mean.mat','d_mean');
sigma_2D = [0.08; 0.02; 0.04; 0.09; 0.01];
% d_real = normrnd(d_mean,sqrt(sigma_2D)*ones(1,Num_t));
% save('d_real.mat','d_real');
load('d_real.mat','d_real');
sigma_D = sqrt(sigma_2D);
Gamma_S = sqrt(Num_j*(1-delta)*(1+Num_j-Num_j*xi)/(1-xi));

%% Parameter for main grid

pi_bf = [25*ones(9,1);35*ones(12,1);25*ones(3,1)]; % dollar/MW
pi_sf = pi_bf-1;

pi_bs = pi_bf;
pi_ss = pi_sf;

E_max = 2.1;

%% Parameter for power system

mpc = loadcase('case5.m');

% Bus
[Nbus,~] = size(mpc.bus);
PD = 0 * mpc.bus(:, 3) / mpc.baseMVA * ones(1, 24); % P demand

% Branch
Ibranch = mpc.branch(:, 1: 2); % branch: from bus, to bus
[Nbranch, ~] = size(Ibranch);
BR_R = mpc.branch(:, 3);
BR_X = mpc.branch(:, 4);
Gbranch = BR_R ./ (BR_R .* BR_R + BR_X .* BR_X);
Bbranch = - BR_X ./ (BR_R .* BR_R + BR_X .* BR_X);
% Sbranch = mpc.branch(:, 6) / mpc.baseMVA; % branch capacity
Sbranch = 10 * ones(size(Gbranch)); %
Sbranch(1) = 4; Sbranch(6) = 2.4;
IFrom = zeros(Nbranch, Nbus);
ITo = zeros(Nbranch, Nbus);
for i = 1: Nbranch
    IFrom(i, Ibranch(i, 1)) = 1;
    ITo(i, Ibranch(i, 2)) = 1;
end

% Generator
Igen = [3; 4; 5];
assert(length(Igen) == Num_i);
Igenbus = zeros(Num_i,Nbus);
for i = 1: Num_i
    Igenbus(i, Igen(i)) = 1;
end

% Prosumer
Ipro = [1; 2; 3; 4; 5];
assert(length(Ipro) == Num_j);
Iprobus = zeros(Num_j,Nbus);
for j = 1: Num_j
    Iprobus(j, Ipro(j)) = 1;
end

%% Parameter for computation

eta_lower = -1000;
M = 1000;
TOL = 0.01;
m_penalty = 1e6;

%% Parameter for linearization

C_m = zeros(Num_j,Num_m);
d_h_m = zeros(Num_j,Num_m);
tau_m = zeros(Num_j,Num_m);
sigma_m = zeros(Num_j,Num_m);
d_e_n = zeros(Num_j,Num_m,Num_t);
for j = 1: Num_j
    tau_m(j,:) = linspace(0.001,0.999,Num_m);
    sigma_m(j,:) = sqrt(sigma_2D(j)./tau_m(j,:)-sigma_2D(j));
    C_m(j,:) = m_C./power(sigma_m(j,:),2);

    d_h_m(j,:) = sqrt(tau_m(j,:).*(1-tau_m(j,:)).*(sigma_2D(j)+power(sigma_m(j,:),2))/(1-delta));
    for t = 1: Num_t
        d_e_n(j,:,t) = (1-tau_m(j,:))*d_mean(j,t)+tau_m(j,:).*(d_real(j,t)+m_sigma*sigma_m(j,:));
    end
end
sigma_mn = zeros(Num_m,Num_m,Num_j,Num_t);
for t = 1: Num_t
    for j = 1: Num_j
        for m = 1: Num_m
            for n = 1: Num_m
                sigma_mn(m,n,j,t) = C_m(j,m)+m_penalty*power(d_e_n(j,n,t)-(1-tau_m(j,m))*d_mean(j,t)-tau_m(j,m)*(d_real(j,t)+m_sigma*sigma_m(j,m)),2);
            end
        end
    end
end
