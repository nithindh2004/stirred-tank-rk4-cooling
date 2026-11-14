function p = default_params()
%DEFAULT_PARAMS  Parameter set identical to src/python/stcool/model.py (SI units).
p.rhoL = 1000; p.V = 0.50; p.ML = p.rhoL*p.V; p.CpL = 4180; p.T0 = 370; p.A = 2.0;
p.Uref = 150; p.Tin = 300; p.tf = 1800;
p.CpC = 4180; p.rhoC = 1000; p.MJ = 10;
p.Nref = 300; p.mref = 0.20; p.eN = 0.25; p.eM = 0.33;
p.Ttarget = 355;
p.D = 0.30; p.Np = 5.0; p.etaAg = 0.75; p.dpRef = 50000; p.etaP = 0.65;
p.rpmMin = 100; p.rpmMax = 600; p.mMin = 0.05; p.mMax = 0.50;
p.flowErrRelSd = 0.05;
end
