%% General parameter
close all
error = 1e-8;

Num_t = 24;

Num_m = 30;

m_C = 0.01;
% m_C = 0.2;
m_sigma = 1; %



delta = 0.95; %
% delta = 0.85;
xi = delta; %

Gamma_T = sqrt(Num_t*(1-delta)*(1+Num_t-Num_t*xi)/(1-xi));

%% Parameter for generator

Num_i = 3;

rho = [30;35;25]*ones(1,Num_t); % dollar/1 MW
rhop = rho;
rhon = [0;0;0]*ones(1,Num_t);
gammap = [0;0;1]*ones(1,Num_t);
gamman = [0;0;1]*ones(1,Num_t);

p_max = [3;2;6]; % MW
% p_min = [0.9;0.6;1.8];
p_min = [1;0.8;2];
r_max = 0.5 * p_max; % reserve
R_max = 0.5 * p_max; % ramp

penalty_shedding = 100; %

%% Parameter for prosumer

Num_j = 5;
d_mean = zeros(5,24);
d_mean(1,:) = (-1)*[0.79,0.73,0.76,0.81,0.90,0.96,0.95,0.89,0.88,0.83,0.87,1.12,...
    1.34,1.52,1.56,1.53,1.39,1.02,0.79,0.71,0.66,0.63,0.60,0.57]-1;
d_mean(2,:) = (-0.5)*[0.59,0.59,0.54,0.51,0.51,0.46,0.41,0.41,0.45,0.45,0.34,0.49,...
    0.79,1.07,1.11,0.99,0.83,0.57,0.35,0.28,0.28,0.33,0.39,0.45]-0.5;
d_mean(3,:) = 2*[0.53,0.55,0.61,0.69,0.71,0.76,0.79,0.83,0.90,0.95,0.96,1.00,...
    0.95,0.90,0.86,0.73,0.71,0.75,0.82,0.95,0.92,0.76,0.67,0.61];
d_mean(4,:) = 2*[0.53,0.55,0.61,0.69,0.71,0.76,0.79,0.83,0.90,0.95,0.96,1.00,...
    0.95,0.90,0.86,0.73,0.71,0.75,0.82,0.95,0.92,0.76,0.67,0.61];
d_mean(5,:) = 3.5*[0.72,0.71,0.68,0.66,0.63,0.61,0.60,0.59,0.60,0.62,0.66,0.80,...
    0.92,1.00,1.04,1.03,0.99,0.86,0.69,0.72,0.62,0.67,0.74,0.82];
figure;
plot(d_mean');
% sigma_2D = [0.04; 0.03; 0.01; 0.04; 0.01];
sigma_2D = [0.08; 0.02; 0.04; 0.09; 0.01];
% sigma_2D = sigma_2D*0.2;
% d_real = normrnd(d_mean,sqrt(sigma_2D)*ones(1,Num_t));
% d_real(1,d_real(1,:)>0) = 0;
% d_real(2,d_real(2,:)>0) = 0;
% save('d_real.mat','d_real');
load('d_real.mat','d_real');
figure;
plot(d_real');
sigma_D = sqrt(sigma_2D);
Gamma_S = sqrt(Num_j*(1-delta)*(1+Num_j-Num_j*xi)/(1-xi));

%% Parameter for power system

mpc = loadcase('case5.m');

% Bus
[Nbus,~] = size(mpc.bus);
PD = [1.5;1.8;0;0;0] * ones(1, 24); % P demand

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
PTDF = zeros(Nbranch+Nbus-1,Nbus);
Arel = IFrom(:,2:Nbus)-ITo(:,2:Nbus);
Mrel = [Arel',zeros(Nbus-1,Nbus-1);diag(BR_X),Arel];
for i = 2:Nbus
    erel = zeros(Nbranch+Nbus-1,1);
    erel(i-1) = 1;
    PTDF(:,i) = Mrel\erel;
end
PTDF = PTDF(1:Nbranch,:);

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
M = 1e6;
TOL = 10;
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
