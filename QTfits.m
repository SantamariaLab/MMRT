function [in] = QTfits(NaCaK,T0, th)
%QTfits returns MMRT and Arrhenius fits of Q10 data and dH using MMRT dCp.
%T0 = 298;      % K
R = 8.31;      % J/(mol*K)
kT = 4.11e-21; 
avagadro = 6.022e23;
inK=@(T) T+273.15;

Q10 = NaCaK(:,2);
Tc = NaCaK(:,1);
T = inK(Tc);

m_fun = fittype('((T+10)./T).*exp((10.*T.*(dH+Cp.*(T-To))+50.*T.*Cp-500.*Cp)./(R.*T.^2.*(T+10)))',...
    'Independent','T','Coefficients',{'Cp','dH'},'Problem',{'To','R'});
a_fun = fittype('exp((10.*G)./(R.*T.*(T+10)))', ...
    'Independent', 'T', 'Coefficients', 'G','Problem',{'R'}); 


m_fitop = fitoptions(m_fun);
a_fitop = fitoptions(a_fun);

m_fitop.Lower=[-1e10 -1e5];
m_fitop.Upper=[1e10 1e5];
m_fitop.StartPoint=[-100 -20];
a_fitop.Lower=[-1e10];
a_fitop.Upper=[1e10];
a_fitop.StartPoint=[-5000];

[sT,sTp]=sort(T);
sQ10=Q10(sTp);

tth = (sT<=th);

[in.m_mdl, in.m_gof] = fit(sT, sQ10, m_fun,'problem',{T0,R}, m_fitop);

[in.a_mdl, in.a_gof] = fit(sT(logical(~tth)), sQ10(logical(~tth)), a_fun,'problem',R, a_fitop);

in.normFitM=norm(in.m_mdl(sT)-sQ10);
in.normFitA=norm(in.a_mdl(sT)-sQ10);


Cp = in.m_mdl.Cp;
T = inK([0:40]);
in.Data=[sT sQ10];
end

