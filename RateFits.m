function [in] = RateFits(rate,T0)
% returns MMRT fits to rate data
% temp=rate(:,1) in C and tau=rate(:,2) in sec
R = 8.31;      % J/(mol*K)
kb = 1.38064852e-23; %m2 kg s-2 K-1
h = 6.62607004e-34; %m2 kg / s
inK=@(T) T+273.15;

T = inK(rate(:,1));
r = 1./rate(:,2); % this variable actually is rise or fall time of AP or time constant for conductances

[sT,sTp]=sort(T);
sr1=r(sTp);

sr=sr1./sr1(1);%normalize to fit the rate coefficient


aafun=fittype('A.*exp(-Ea./(R.*T))',...
    'Independent','T','Coefficients',{'A','Ea'},'Problem',{'R'});
mmfun=fittype('(kb.*T/h).*exp(-(dH+Cp.*(T-To))./(R.*T)+(dS + Cp.*(log(T)-log(To)))./R)',...
    'Independent','T','Coefficients',{'Cp','dH','dS'},'Problem',{'To','R','kb','h'});

aafitop = fitoptions(aafun);
mmfitop = fitoptions(mmfun);
% mmfitop.StartPoint=[-2050 -600 1120];
%mmfitop.Lower=[-inf -inf 0]
%mmfitop.Upper=[0 0 inf]

[in.aa_mdl, in.aa_gof] = fit(sT,sr,aafun,'problem',{R},aafitop);
[in.mm_mdl, in.mm_gof] = fit(sT,sr,mmfun,'problem',{T0,R,kb,h},mmfitop);

in.normFitA=norm(in.aa_mdl(T)-r);
in.normFitM=norm(in.mm_mdl(T)-r); 
in.Data=[sT sr];
in.DataOrig=[sT sr1];
end

