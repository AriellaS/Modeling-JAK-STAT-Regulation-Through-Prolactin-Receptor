%% Run model simulations
    % Input: parameter values and initial conditions
    % Output: predicted time courses for all molecular species

clear all
close all

load fittedParams;

% Parameter values
params(2,1) = 0.000056;
params(3,1) = 0.0056;
params(4,1) = fitted_params.k3;
params(5,1) = 0.2;
params(6,1) = 0.005;
params(7,1) = fitted_params.k5;
params(8,1) = 0.8;
params(9,1) = 0.4;
params(10,1) = 0.000256721177985165;
params(11,1) = fitted_params.deg_ratio;
params(12,1) = fitted_params.k8A;
params(13,1) = 0.1;
params(14,1) = fitted_params.k8A * fitted_params.mult8B;
params(15,1) = 0.1;
params(16,1) = fitted_params.k8A * fitted_params.mult8AB;
params(17,1) = 0.1;
params(18,1) = fitted_params.k9;
params(19,1) = 0.2;
params(20,1) = 0.003;
params(21,1) = fitted_params.k11;
params(22,1) = 0.2;
params(23,1) = 0.003;
params(24,1) = 0.0000002;
params(25,1) = fitted_params.k_13;
params(26,1) = fitted_params.k14A;
params(27,1) = fitted_params.k14A * fitted_params.mult14B;
params(28,1) = fitted_params.k14A * fitted_params.mult14AB;
params(29,1) = fitted_params.k15;
params(30,1) = 0.2;
params(31,1) = fitted_params.k16;
params(32,1) = 0.0355;
params(33,1) = fitted_params.k17outA;
params(34,1) = 0.0355;
params(35,1) = fitted_params.k17outA * fitted_params.mult17B;
params(36,1) = 0.01;
params(37,1) = 400;
params(38,1) = fitted_params.k19;
params(39,1) = 0.01;
params(40,1) = fitted_params.k21;
params(41,1) = 0.1;
params(42,1) = fitted_params.k22;
params(43,1) = fitted_params.k23;
params(44,1) = fitted_params.k24;
params(45,1) = fitted_params.k25a;
params(46,1) = 400;
params(47,1) = 0.001;
params(48,1) = fitted_params.k27;
params(49,1) = fitted_params.k28;
params(50,1) = 0.01;
params(51,1) = 0.5;
params(52,1) = 1.2;
params(53,1) = 1.36;
params(54,1) = fitted_params.totalSTAT;
params(55,1) = fitted_params.k30a;
params(56,1) = 400;
params(57,1) = 0.001;
params(58,1) = 0.0005;
params(59,1) = 0.01;
params(60,1) = 1.92540883488874E-05;
params(61,1) = 0;

% Initial values
initvalues(1,1)= 9.09;
initvalues(2,1)= fitted_params.RJ;
initvalues(3,1)= 0;
initvalues(4,1)= 0;
initvalues(5,1)= fitted_params.SHP2;
initvalues(6,1)= fitted_params.PPX;
initvalues(7,1)= fitted_params.PPN;
initvalues(8,1)= 0;
initvalues(9,1)= 0;
initvalues(10,1)= 0;
initvalues(11,1)= 0;
initvalues(12,1)= 0;
initvalues(13,1)= 0;
initvalues(14,1)= 0;
initvalues(15,1)= 0;
initvalues(16,1)= 0;
initvalues(17,1)= 0;
initvalues(18,1)= 0;
initvalues(19,1)= 0;
initvalues(20,1)= 0;
initvalues(21,1)= 0;
initvalues(22,1)= 0;
initvalues(23,1)= 0;
initvalues(24,1)= 0;
initvalues(25,1)= 0;
initvalues(26,1)= 0;
initvalues(27,1)= 0;
initvalues(28,1)= 0;
initvalues(29,1)= 0;
initvalues(30,1)= 0;
initvalues(31,1)= 0;
initvalues(32,1)= 0;
initvalues(33,1)= 0;
initvalues(34,1)= 0;
initvalues(35,1)= 0;
initvalues(36,1)= 0;
initvalues(37,1)= 0;
initvalues(38,1)= 0;
initvalues(39,1)= 0;
initvalues(40,1)= 0;
initvalues(41,1)= 0;
initvalues(42,1)= 0;
initvalues(43,1)= 0;
initvalues(44,1)= 0;
initvalues(45,1)= 0;
initvalues(46,1)= 0;
initvalues(47,1)= 0;
initvalues(48,1)= 0;
initvalues(49,1)= 0;
initvalues(50,1)= 0;
initvalues(51,1)= 0;
initvalues(52,1)= 0;
initvalues(53,1)= 0;
initvalues(54,1) = 0;
initvalues(55,1) = 0;
initvalues(56,1) = 50;

% Set number of samples
n_samples = 1E4;

% Set options
options = odeset('RelTol',1e-9,'AbsTol',1e-12,'NonNegative',[1:length(initvalues)]);

% Set runtime
predTime = [0:60:24*3600];

% creates a matrix where each row is a sample and each column is a nonzero init value
samples = zeros(n_samples, 6);
for n = 1 : n_samples
	samples(n,1) = random(initvalues(1,1));
	samples(n,2) = random(initvalues(2,1));
	samples(n,3) = random(initvalues(5,1));
	samples(n,4) = random(initvalues(6,1));
	samples(n,5) = random(initvalues(7,1));
	samples(n,6) = random(initvalues(56,1));
end

% simulate the model
parpool(4);
results = zeros(n_samples, length(predTime), 61);
parfor n = 1 : n_samples
	disp(n); % just to see how fast the simulation is running
	results(n,:,:) = simulator(predTime,options,params,initvalues,samples,n);
end
