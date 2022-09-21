function [in] = logRateFits(rate,T0)
% returns MMRT fits to rate data
% temp=rate(:,1) in C and tau=rate(:,2) in sec
R = 8.31;      % J/(mol*K)
kb = 1.38064852e-23; %m2 kg s-2 K-1
h = 6.62607004e-34; %m2 kg / s
inK=@(T) T+273.15;

T = inK(rate(:,1));
% rate(:,1)
% [Tnorp]=find((rate(:,1)>20).*(rate(:,1)<26));
% if isempty(Tnorp)
%     [mv,Tnorp]=max(rate(:,1))
% end
[mv,Tnorp]=min(rate(:,1))
r = 1./rate(:,2);
r=r./r(Tnorp(1));%normalize to about the same temperature all datasets


[sT,sTp]=sort(T);
sr=r(sTp);

srn=sr./sr(1);%normalize to fit the rate coefficient

srnl=log(srn);

aafun=fittype('A - Ea./(R.*T)',...
    'Independent','T','Coefficients',{'A','Ea'},'Problem',{'R'});

% mmfun=fittype('(kb.*T/h).*exp(-(dH+Cp.*(T-To))./(R.*T)+(dS + Cp.*(log(T)-log(To)))./R)',...
%     'Independent','T','Coefficients',{'Cp','dH','dS'},'Problem',{'To','R','kb','h'});
mmfun=fittype('log(kb.*T/h)+(-(dH+Cp.*(T-To))./(R.*T)+(dS + Cp.*(log(T)-log(To)))./R)',...
    'Independent','T','Coefficients',{'Cp','dH','dS'},'Problem',{'To','R','kb','h'});

% mmfun=fittype('(kb.*T/h).*exp(-(Cp.*(T-T0))./(R.*T)-dG./(R.*T)+Cp.*log(T./T0)./R)',...
%     'Independent','T','Coefficients',{'Cp','dG'},'Problem',{'T0','R','kb','h'});
%mmfun=fittype('log(kb.*T/h)-(Cp.*(T-T0))./(R.*T)-dG./(R.*T)+(Cp.*(log(T)-log(T0)))./R',...
%    'Independent','T','Coefficients',{'Cp','dG'},'Problem',{'T0','R','kb','h'});
% mmfun=fittype('log(kb./h)+log(T).*(1+Cp./R)+Cp.*(-1/R+T0/(R.*T))-dG/R+(Cp.*(-log(T0)))./R',...
%     'Independent','T','Coefficients',{'Cp','dG'},'Problem',{'T0','R','kb','h'});
aafitop = fitoptions(aafun);
mmfitop = fitoptions(mmfun);
mmfitop.StartPoint = [-1375 100 2];
% mmfitop.Weights=ones(1,length(r));
% mmfitop.Weights(1)=100;
% mmfitop.Weights(2:end-1)=20;
% mmfitop.Weights(end)=100;
%mmfitop.Lower =[-Inf 0 -Inf];
%mmfitop.Upper =[0 Inf Inf];

% mmfitop.StartPoint=[-2050 -600 1120];
%mmfitop.Lower=[-inf -inf 0]
%mmfitop.Upper=[0 0 inf]

[in.aa_mdl, in.aa_gof] = fit(sT,srnl,aafun,'problem',{R},aafitop);
[in.mm_mdl, in.mm_gof] = fit(sT,srnl,mmfun,'problem',{T0,R,kb,h},mmfitop);
[in.lin_mdl,in.lin_gof]= fit(sT,srnl,'poly1');
% clf
% plot(T,r);
% hold on
% plot(T,in.mm_mdl(T),'k')
in.normFitA=norm(in.aa_mdl(T)-r);
in.normFitM=norm(in.mm_mdl(T)-r); 
in.Data=[sT sr];
end

