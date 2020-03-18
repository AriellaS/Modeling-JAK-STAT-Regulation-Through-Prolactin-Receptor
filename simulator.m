function result = simulator(predTime,options,params,initvalues,samples,n)

	% assign free variables back to initvalues
	initvalues(1,1) = samples(n,1);
	initvalues(2,1) = samples(n,2);
	initvalues(5,1) = samples(n,3);
	initvalues(6,1) = samples(n,4);
	initvalues(7,1) = samples(n,5);
	initvalues(56,1) = samples(n,6);

	% Assign variable names for calculated quantities
	kdeg = params(10,1);
	RJ = initvalues(2,1);
	Vratio = params(51,1);
	ncratioA = params(52,1);
	ncratioB = params(53,1);
	totalSTAT = params(54,1);
	k17outA = params(33,1);
	k17outB = params(35,1);
	k34 = params(60,1);
	BCL = initvalues(56,1);

	% Calculate quantities
	ksyn = kdeg*RJ; % syn rate of RJ
	S5Ac = totalSTAT/(1 + ncratioA*Vratio); % STAT5 in cytosol
	S5An = (totalSTAT - S5Ac)/Vratio; %STAT5 in nucleus
	S5Bc = totalSTAT/(1 + ncratioB*Vratio); % STAT5 in cytosol
	S5Bn = (totalSTAT - S5Bc)/Vratio; %STAT5 in nucleus
	k17inA = S5An/S5Ac*(Vratio)*k17outA;
	k17inB = S5Bn/S5Bc*(Vratio)*k17outB;
	k35 = k34*BCL;

	%Assign them back to params and initvalues
	params(1,1) = ksyn;
	initvalues(3,1) = S5Ac;
	initvalues(35,1) = S5An;
	initvalues(4,1) = S5Bc;
	initvalues(36,1) = S5Bn;
	params(32,1) = k17inA;
	params(34,1) = k17inB;
	params(61,1) = k35;

	% run odes
	[~, predConc] = ode15s(@core_file,predTime,initvalues,options,params);

	% Calculated quantities
	total_pStatA = predConc(:,13) +2.*predConc(:,15) + predConc(:,16) + predConc(:,19) + 2*predConc(:,21) + predConc(:,23) + predConc(:,24) + predConc(:,27) + 2*predConc(:,28) + predConc(:,29) + predConc(:,31) + predConc(:,33) + 2*predConc(:,37) + predConc(:,38) + predConc(:,40) + predConc(:,43);
	total_pStatB = predConc(:,14) +2.*predConc(:,17) + predConc(:,16) + predConc(:,20) + 2*predConc(:,22) + predConc(:,23) + predConc(:,25) + predConc(:,26) + 2*predConc(:,30) + predConc(:,29) + predConc(:,32) + predConc(:,34) + 2*predConc(:,39) + predConc(:,38) + predConc(:,41) + predConc(:,42);

	pStatA_norm = total_pStatA./total_pStatA(31);% normalized to 30 minute
	pStatB_norm = total_pStatA./total_pStatA(31);% normalized to 30 minute

	Stat_cytoA = predConc(:,3) + predConc(:,11) + predConc(:,13) + 2*predConc(:,15) + predConc(:,16) + predConc(:,19) + 2*predConc(:,21) + predConc(:,23) + 2*predConc(:,24) + predConc(:,26) + predConc(:,27) + predConc(:,48);
	Stat_cytoB = predConc(:,4) + predConc(:,12) + predConc(:,14) + 2*predConc(:,17) + predConc(:,16) + predConc(:,20) + 2*predConc(:,22) + predConc(:,23) + 2*predConc(:,25) + predConc(:,26) + predConc(:,27) + predConc(:,49);

	Stat_nucleusA = 2*predConc(:,28) + predConc(:,29) + predConc(:,31) + predConc(:,33) + predConc(:,35) + 2*predConc(:,37) + predConc(:,38) + 2*predConc(:,40) + predConc(:,42) + predConc(:,43);
	Stat_nucleusB = 2*predConc(:,30) + predConc(:,29) + predConc(:,32) + predConc(:,34) + predConc(:,36) + 2*predConc(:,39) + predConc(:,38) + 2*predConc(:,41) + predConc(:,42) + predConc(:,43);

	nucleus_cyto_ratioA = Stat_nucleusA./Stat_cytoA;
	nucleus_cyto_ratioB = Stat_nucleusB./Stat_cytoB;

	% Bcl-xL fold change
	Bcl = predConc(:,56)./predConc(1,56);

	% Create results matrix
	result = predConc;
	result(:,57) = total_pStatA;
	result(:,58) = total_pStatB;
	result(:,59) = nucleus_cyto_ratioA;
	result(:,60) = nucleus_cyto_ratioB;
	result(:,61) = Bcl;

end
