%% Figure 1A,B,C,D
clear
load Q10z_Fidel

inK = @(T) T+273.15;


T0=25;
INa=QTfits(Q10.Na,inK(T0),inK(0)); 
IK=QTfits(Q10.K,inK(T0),inK(0));
ICa=QTfits(Q10.Ca,inK(T0),inK(0));

figure(1)
clf
subplot 311
plot((INa.Data(:,1))-273,INa.Data(:,2),'o','MarkerSize',3,'MarkerEdgeColor',[0 0 0])
hold on
plot([5:40],INa.m_mdl(inK([5:40])),'k')
plot([5:40],INa.a_mdl(inK([5:40])),'k--')
box off
legend('exp','MMRT','Arrhenius') 
legend boxoff
ylabel('Q10')
title('Na+')
ylim([0 10])

subplot 312
plot((IK.Data(:,1)-273),IK.Data(:,2),'o','MarkerSize',3,...
    'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
hold on
plot([5:40],IK.m_mdl(inK([5:40])),'k')
plot([5:40],IK.a_mdl(inK([5:40])),'k--')
box off
title('K+')
ylabel('Q10')
ylim([0 10])

subplot 313
plot((ICa.Data(:,1)-273), ICa.Data(:,2), 'sq','MarkerSize',3,'MarkerEdgeColor',[0 0 0])
hold on 
plot([5:40],ICa.m_mdl(inK([5:40])),'k')
plot([5:40],ICa.a_mdl(inK([5:40])),'k--')
box off
xlabel('T (C)')
title('Ca++')
ylabel('Q10')
ylim([0 10])

%% Figure 1D
%% Rate fit MMRT
clear
load Q10z_Bahram r

fitingdata=["r.Ca.Ca1" "r.Ca.Ca4" "r.Ca.Ca5" "r.Ca.Ca51" "r.Ca.Ca6" "r.K.K1" "r.K.K11" "r.K.K21" "r.K.K23" "r.K.K22" ...
    "r.Na.Na31" "r.Na.Na61" "r.Na.Na62" "r.Na.Na71" "r.Na.Na82"];


inK=@(T) T+273.15;
T0=inK(25);
T2=100;

% % log fit
clf
figure(10)
clear logfitsR
c=1;
for i=[1:15] 
    d=eval(fitingdata(i));
    dummy=logRateFits(d,T0);
    if i<=5
        markerS='o';
        displayN='Ca';
        markerF=[1 1 1];
        linS='k';
    elseif i<=10
        markerS='o';
        markerF=[0 0 0];
        displayN='K';
        linS='k'; %'k--';
    else
        displayN='Na';
        markerF=[1 1 1 ];
        markerS='sq';
        linS='k'; %'k-.';
    end
    semilogy((dummy.Data(:,1))-273,(dummy.Data(:,2)),markerS,'markersize',3,...
            'MarkerFaceColor',markerF,'MarkerEdgeColor',[0 0 0])
        hold on
    semilogy([0:T2],exp(dummy.mm_mdl(inK([0:T2]))).*dummy.Data(1,2),...
            linS,'DisplayName',displayN)
    logfitsR(c)=dummy;
    c=c+1;
%     pause
end
box off
xlabel('T (C)')
ylabel('rate (1/s)')


%Creating table
R = 8.31;      % J/(mol*K)
cint95=[];ratetable=[];Cps=[];dHs=[];dSs=[];rsqs=[];Top=[];
for i=1:length(logfitsR)
    Cps(i)=logfitsR(i).mm_mdl.Cp;
    dHs(i)=logfitsR(i).mm_mdl.dH;
    dSs(i)=logfitsR(i).mm_mdl.dS;
    rsqs(i)=logfitsR(i).mm_gof.rsquare;
    if length(logfitsR(i).Data)<=3
        dummy=[0 0 0; 0 0 0];
    else
        dummy=confint(logfitsR(i).mm_mdl)-coeffvalues(logfitsR(i).mm_mdl);
    end
    cint95{i}=dummy(2,:);
    cps95(i)=dummy(2,1);
    dhs95(i)=dummy(2,2);
    dss95(i)=dummy(2,3);
    Top(i)=(Cps(i).*inK(25)-dHs(i))./(Cps(i)+R);
    datapoints(i)=size(logfitsR(i).Data,1);
    diff(sign(diff(logfitsR(i).mm_mdl(inK(0:5000)))));
    ctype(i)=extractBetween(fitingdata{i},'.','.')
end
ratetable=table([1:15]',datapoints',ctype',Cps',cps95',dHs',dhs95',dSs',dss95',rsqs',Top')
ratetable.Properties.VariableNames={'OrigN','DP','Cond','Cps','cps95','dHs','dhs95','dSs','dSs95','rsqs','Top'}

bodyT=[37,38,37,37,39,37,37,37,37,37,21,21,21,37,7];
bodyT2=bodyT([1 3:5 8:13 14 15]); 
%If a line fits better the data then we discard them. 
rsqcomp=[]
for i=1:15
    dummy=logfitsR(i);
    rounv=round(100.*[dummy.mm_gof.rsquare dummy.lin_gof.rsquare])./100;
    drv=100*(diff(rounv)./rounv(2));
    rsqcomp(i,:)=[i  rounv drv dummy.mm_gof.adjrsquare dummy.lin_gof.adjrsquare]
end

tablefil1=ratetable(logical(rsqcomp(:,4)),:);
tablefil2=tablefil1([1:12],:);
NaST=tablefil2(logical(strcmp(tablefil2.Cond,'Na')),:);
CaST=tablefil2(logical(strcmp(tablefil2.Cond,'Ca')),:);
KST=tablefil2(logical(strcmp(tablefil2.Cond,'K')),:);

figure(2)
subplot(1,2,1)
cla;
T2=80;
for a=1:length(tablefil2.Top)
    dummy=logfitsR(tablefil2.OrigN(a));
    condT=tablefil2.Cond{a};
    if strcmp(condT,'Ca')
        markerS='sq';
        displayN='Ca';
        markerF=[1 1 1];
        linS='k';
    elseif strcmp(condT,'K')
        markerS='o';
        markerF=[0 0 0];
        displayN='K';
        linS='k'; %'k--';
    else
        displayN='Na';
        markerF=[1 1 1 ];
        markerS='o';
        linS='k'; %'k-.';
    end
    semilogy((dummy.Data(:,1))-273,(dummy.Data(:,2)),markerS,'markersize',7,...
            'MarkerFaceColor',markerF,'MarkerEdgeColor',[0 0 0])
        hold on
    semilogy([0:T2],exp(dummy.mm_mdl(inK([0:T2]))).*dummy.Data(1,2),...
            linS,'DisplayName',displayN)
end
box off
xlim([0 50])
ylim([1 100])
ylabel('Rate coefficient (k(T))')
xlabel('Temperature (^{o}C)')

set(2,'WindowStyle','normal')
exportgraphics(gcf,'MMRT-rateCoeff.eps')
set(2,'WindowStyle','docked')

%meanTop=[mean(Top([1 3:5])) mean(Top(8:10)) mean(Top([11:13 16]))]-273

%% Figure 2
subplot(2,2,3)
cla
q10C=1e3.*[-3.61 	94.78;
-4.15 	86.5;
-2.49 	76.72 ];  %Na,K,Ca
Topq10=(q10C(:,1).*inK(25)-q10C(:,2))./(q10C(:,1)+R)-273

TopNa=NaST.Top;
[TopsNa,ToppNa]=sort(TopNa);

TopK=KST.Top;
[TopsK,ToppK]=sort(TopK);

TopCa=CaST.Top;
[TopsCa,ToppCa]=sort(TopCa);

% cla
% h=plot([ [TopsNa' TopsK' TopsCa']-273],'o'); %without number 6
% h1=h.Parent;
% h1.XTick=[1:length([TopNa' TopK' TopCa'] )];
% h1.XTickLabel={NaST.Cond{ToppNa},...
%     KST.Cond{ToppK},CaST.Cond{ToppCa}};%,fitingdata(Topp)}
% h1.XTickLabelRotation=80;
% hold on
% box off
% plot(3,Topq10(1),'ksq','MarkerSize',10);
% plot(8,Topq10(2),'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10)
% plot(10.5,Topq10(3),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',10)

cla
h=errorbar([1 2 3],[mean(TopNa) mean(TopK) mean(TopCa)]-273,...
    [std(TopNa)/sqrt(5) std(TopK)./sqrt(3) std(TopCa)./sqrt(4)],'k.') %std(TopNa([1:4 6]))/sqrt(5) to get rid of 250 C
h1=h.Parent;
h1.XTick=[1 2 3];
h1.XTickLabel={'Na','K','Ca'}
xlim([0.5 3.5])
box off
hold on
% scatter([1 1 1 1 1 1 2 2 2 3 3 3 3],[TopNa' TopK' TopCa']-273,'k')
scatter([1 1 1 1 1 ],[TopNa' ]-273,'ks')
scatter([2 2 2 ],[ TopK']-273,'k','filled')
scatter([ 3 3 3 3],[ TopCa']-273,'k')
plot(1.2,Topq10(1),'ksq','MarkerSize',10);
plot(2.2,Topq10(2),'ko','MarkerFaceColor',[0 0 0],'MarkerSize',10)
plot(3.2,Topq10(3),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',10)
%'NaQ10','KQ10','CaQ10' 
tablefil3=tablefil2(logical(tablefil2.Top<(150+273)),:);
writetable(tablefil3,'mmrt-rate.xlsx') %final table 
ylabel('T_{opt}')


[mean(TopNa)-273 std(TopNa)./sqrt(numel(TopNa))]
[mean(TopK)-273 std(TopK)./sqrt(numel(TopK))]
 [mean(TopCa)-273 std(TopCa)./sqrt(numel(TopCa))]





%% Figure 3
R = 8.31;      % J/(mol*K)
kb = 1.38064852e-23; %m2 kg s-2 K-1
hplank = 6.62607004e-34; %m2 kg / s
inK=@(T) T+273.15;
T0=inK(25);
T2=100;
T=inK(5:T2);

%from MMRT Q10 from above these are the coefficients for easy use
CpNa=-2860; dHNa=9.011e+04;
CpK=-4144; dHK=8.898e+04;
CpCa=-3607; dHCa=9.718e+04;

%Based on tref k(tref) = 1, so we adjusted dS in each MMRT-k(T). 
tref=20;
dSNa=63.1; dSK=59.5; dSCa=87.4;

mmfunK=@(T) (kb.*T/hplank).*exp(-(dHK+CpK.*(T-T0))./(R.*T)+(dSK+CpK.*(log(T)-log(T0)))./R);
mmfunNa=@(T) (kb.*T/hplank).*exp(-(dHNa+CpNa.*(T-T0))./(R.*T)+(dSNa+CpNa.*(log(T)-log(T0)))./R);
mmfunCa=@(T) (kb.*T/hplank).*exp(-(dHCa+CpCa.*(T-T0))./(R.*T)+(dSCa+CpCa.*(log(T)-log(T0)))./R);

% % % plotting Na,K,Ca
figure(3)
clf
plot(T-273.15,mmfunNa(T)./mmfunNa(inK(tref)),'k')
hold on
plot(T-273.15,mmfunK(T)./mmfunK(inK(tref)),'k')
plot(T-273.15,mmfunCa(T)./mmfunCa(inK(tref)),'k')
box off

xlabel('T (C)')
ylabel('rate coefficient')

%% Figure 4 _ HH spikes
clear
tstop  = 100; %ms
tstep  = 1e-3; %ms
mVinit = -65; %mV

Ncells = 1; %only modeling one cell
Cm=1;  %microF/cm^2
dt =tstep;
t=0:dt:tstop; %ms
v0 = mVinit*ones([1 Ncells]);              % mV  initial value
vrest=v0(1,1);% mV  the resting  potential.

gK=36; gNa=120; gL=0.3;      % channel conductances: mS/cm2 -- (micro A/mV)/cm^2
EK=-12 + v0(1,1); ENa=115 + v0(1,1); EL=10.6 + v0(1,1);  % channel reversal potentials: mV
m=0.0529;h=0.596; n=0.3177; % Initial vlaues:  steady state for 0 input n=0.3177, m=0.0529, h=0.596, v=-64.9997

NetProp.Ncells=Ncells;

NetProp.dt=dt;
NetProp.Cm=Cm;
NetProp.v0=v0;
NetProp.vrest=vrest;
NetProp.gK=gK;
NetProp.gNa=gNa;
NetProp.gL=gL;
NetProp.EK=EK;
NetProp.ENa=ENa;
NetProp.EL=EL;
NetProp.m=m;
NetProp.h=h;
NetProp.n=n;
NetProp.Noise=0;

NetProp.MemoryWindow=1000/dt;

Iamp_1=[1:12];
Iamp_2=[13:24];
Iamp=[Iamp_1 Iamp_2];

I=ones(length(t),Ncells);
%I(1:10/dt,Ncells)=0;
f2safe='HH_fracN_firingrateVI';
IampV=[1:1:2000];
IampVs=[11.*ones(1,6) 13.*ones(1,5) 15.*ones(1,5) 18.*ones(1,5) 20.*ones(1,5) 27.*ones(1,5)]% 24.*ones(1,5)]; %this accelerates the simulations because we already did the hard work of finding the rheobase
 
Tvector=[0:30];
Vth=-20; %-65+45
for model="MMRT-HK3"%the other address different models and approaches we took ["Q10_3","ARRHENIUS","MMRT","MMRT-K","MMRT-HK","MMRT-HK2",]
    for a=1:length(Tvector)
        foundspike=0;
        c=1;
        IampV=[IampVs(a):2000];
        disp(['Temp: ' num2str(Tvector(a))])
        while ~foundspike
            disp(['      Iamp: ' num2str(IampV(c))]);
            I2inj=0*ones(size(t));
            I2inj(40/dt:40.5/dt)=IampV(c);
            out02(1)=HHclassic(NetProp,I2inj,t,Tvector(a),model);
            plot(out02.v)
            foundspike=(sum(out02.v>Vth)>0);
            lastI=IampV(c);
            c=c+1;
        end
        if model=='Q10_3'
            spikeVsI3(a).st=out02; spikeVsI3(a).T=Tvector(a); spikeVsI3(a).Iamp=lastI;
        elseif model=='ARRHENIUS'
            spikeVsIA(a).st=out02; spikeVsIA(a).T=Tvector(a); spikeVsIA(a).Iamp=lastI;
        elseif model=='MMRT'
            spikeVsIM(a).st=out02; spikeVsIM(a).T=Tvector(a); spikeVsIM(a).Iamp=lastI;
        elseif model=='MMRT-K'
            spikeVsIMK(a).st=out02; spikeVsIMK(a).T=Tvector(a); spikeVsIMK(a).Iamp=lastI;
        elseif model=='MMRT-HK'
            spikeVsIMHK(a).st=out02; spikeVsIMHK(a).T=Tvector(a); spikeVsIMHK(a).Iamp=lastI;
        elseif model=='MMRT-HK2'
            spikeVsIMHK2(a).st=out02; spikeVsIMHK2(a).T=Tvector(a); spikeVsIMHK2(a).Iamp=lastI;
        elseif model=='MMRT-HK3'
            spikeVsIMHK3(a).st=out02; spikeVsIMHK3(a).T=Tvector(a); spikeVsIMHK3(a).Iamp=lastI;
        end
    end
end
%Uncomment this to save the data
save HHsimulfinal
% save('HHsimul.mat','-v7.3')
% save('HHsimulMHK.mat','-v7.3')
% save('HHsimulMHK2.mat','-v7.3')

%% Figure 4A
% this is to
clear
inK=@(T) T+273.15;
Tvector=[0:35];
T0=inK(25);
T2=100;

%Data stracted from Hodgkin-Katz (fall and rise times in ms)
HKfall=[4.518117483 5.893240005 6.596273735 7.603470868 9.388035654 11.45826295 12.96035883 17.24419531 17.71024148 18.14005675 18.60478142 20.04234444 27.79300499 29.93943833 32.53473844;...
     3.116689052 2.498924633 2.223011922 1.781210109 1.338323119 0.992939111 0.855986365 0.537485403 0.503728903 0.491041686 0.463234374 0.415350062 0.2126709 0.189679315 0.172679559];
HKrise=[3.241970157 4.488252342 4.675946113 5.982121091 6.754604529 7.725468927 9.475496256 11.52688007 13.11782413 17.26256125 17.679139 18.10696056 18.54680522 20.00404808 27.6852431 29.9016104 32.47599857;...
     1.303773066 1.206589723 1.114534291 1.017955873 1.019365913 0.894125655 0.759716548 0.690212275 0.634927534 0.490428232 0.453197385 0.456567175 0.433283491 0.380381245 0.271222615 0.254804536 0.222658422];

HKriserate=RateFits(HKrise',T0);
HKfallrate=RateFits(HKfall',T0);

figure(4)
clf
t=tiledlayout(5,2)
%Plot original HK rates and fits. Please pay attention that the fit is to
%k(T) so to replot the r(T)=k(T).*k(tref), which was the lowest temperatuer
%in each dataset
nexttile([1 1])
plot(HKriserate.Data(:,1)-273.15,HKriserate.Data(:,2).*HKriserate.DataOrig(1,2),'o','MarkerSize',5,'MarkerEdgeColor',[0 0 0])
hold on
plot(HKfallrate.Data(:,1)-273.15,HKfallrate.Data(:,2).*HKfallrate.DataOrig(1,2),'o','MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
plot(HKriserate.mm_mdl(inK(Tvector)).*HKriserate.DataOrig(1,2),'k')
plot(HKfallrate.mm_mdl(inK(Tvector)).*HKfallrate.DataOrig(1,2),'k--')
box off
legend('rise','fall','H&K rise rate','H&K fall rate')
legend boxoff
ylabel('rise and fall rate (1/ms)')
xlabel('T (C)')

%% Figure 4B
% first you have save the data from above
 load HHsimulfinal
t=0:1e-3:100;
s=241:248;
T=0:5:35;

nexttile

for i=[1 3 5 7] %2:8 %plots more examples that in the Figure
    %subplot(s(i))
    nexttile
    title(sprintf('%i C', T(i)))
    hold on
    %plot(t,spikeVsI3(5*(i-1)+1).st.v,'k :')
    %plot(t,spikeVsIA(5*(i-1)+1).st.v,'k--')
    plot(t-35,spikeVsIMHK3(5*(i-1)+1).st.v,'k')

    if mod(i,4)==1, ylabel('V (mV)'), end
    if i>4, xlabel('t (ms)'), end
    if i==4, legend({'Q10 = 3', 'Arrhenius', 'MMRT', 'MMRT-K', 'MMRT-HK', 'MMRT-HK2'}), legend boxoff, end
    axis([-5 30 -80 50])
end
%% Figure 4 C
%rise/fall time
% first we analize the MMRT-HH action potentials
%plot in a different figure and come back
figure(10)
Tvector=[0:30];
t=0:1e-3:100;
%Analize rise/fall times of the simulated MMRT-HH spikes
for i=1:length(Tvector)
    spkMHK3=[t',spikeVsIMHK3(i).st.v];
    [spikeVsIMHK3(i).rise,spikeVsIMHK3(i).fall]=apdur(spkMHK3(35000:70000,:),i);
    riseMHK3(i)=spikeVsIMHK3(i).rise;
    fallMHK3(i)=spikeVsIMHK3(i).fall;
end

figure(4)


nexttile([1 2])
%this are the same fits for 4A but 1/rate to plote the fall and rise times
plot(1./(HKriserate.mm_mdl(inK(Tvector)).*HKriserate.DataOrig(1,2)),'k')
hold on
plot(1./(HKfallrate.mm_mdl(inK(Tvector)).*HKfallrate.DataOrig(1,2)),'k--')

plot(Tvector,riseMHK3,'k-') %rise measured from MMRT-HH action potentials
plot(Tvector,fallMHK3,'k-') %rise measured from MMRT-HH action potentials

%Hodgkin-Katz
plot([4.518117483 5.893240005 6.596273735 7.603470868 9.388035654 11.45826295 12.96035883 17.24419531 17.71024148 18.14005675 18.60478142 20.04234444 27.79300499 29.93943833 32.53473844],...
     [3.116689052 2.498924633 2.223011922 1.781210109 1.338323119 0.992939111 0.855986365 0.537485403 0.503728903 0.491041686 0.463234374 0.415350062 0.2126709 0.189679315 0.172679559],...
     'o','MarkerSize',3,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]) %'og')%fall
plot([3.241970157 4.488252342 4.675946113 5.982121091 6.754604529 7.725468927 9.475496256 11.52688007 13.11782413 17.26256125 17.679139 18.10696056 18.54680522 20.00404808 27.6852431 29.9016104 32.47599857],...
     [1.303773066 1.206589723 1.114534291 1.017955873 1.019365913 0.894125655 0.759716548 0.690212275 0.634927534 0.490428232 0.453197385 0.456567175 0.433283491 0.380381245 0.271222615 0.254804536 0.222658422],...
     'o','MarkerSize',3,'MarkerEdgeColor',[0 0 0])%'*g')%rise
box off
ylabel('rise and fall time (ms)')
xlabel('T (C)')


%% Figure 4 D
%Q10 from the fits to rise and fall of the HK data
HKrQ10=HKriserate.mm_mdl(inK(Tvector(11:31)))./HKriserate.mm_mdl(inK(Tvector(1:21))); 
HKfQ10=HKfallrate.mm_mdl(inK(Tvector(11:31)))./HKfallrate.mm_mdl(inK(Tvector(1:21))); 
% Calculate Q10 based on simulation of MMRT-HH action potentials
for i=1:length(Tvector)-10
    spikeVsIMHK3(i).rQ10=spikeVsIMHK3(i).rise/spikeVsIMHK3(i+10).rise;
    spikeVsIMHK3(i).fQ10=spikeVsIMHK3(i).fall/spikeVsIMHK3(i+10).fall;
    q10rMHK3(i)=spikeVsIMHK3(i).rQ10;
    q10fMHK3(i)=spikeVsIMHK3(i).fQ10;
end


nexttile
plot(HKrQ10,'k:'); %based on fits to HK data 
hold on
plot(Tvector(1:21),q10rMHK3,'k-')
plot([5 10 20],[2.7 2.07 1.54],'ok','MarkerSize',3,'MarkerEdgeColor',[0 0 0])
plot([0 10 20],[2.1 1.78 1.42],'^','MarkerSize',3,'MarkerEdgeColor',[0 0 0])
plot([0 10 20],[2.56 1.85 1.54],'sq','MarkerSize',3,'MarkerEdgeColor',[0 0 0])
plot([0 10 20],[2.66 1.85 1.63],'+','MarkerSize',3,'MarkerEdgeColor',[0 0 0])
plot([0 10 20],[3.16 1.78 1.54],'d','MarkerSize',3,'MarkerEdgeColor',[0 0 0]);
box off
ylim([0 9])
xlabel('T ref (C)')
ylabel('Q10 rise')
legend({'MMRT', 'HH', 'L opal', 'L. pealei', 'L. plei', 'S. sepioidea','HK fit'})
legend boxoff 

nexttile
plot(HKfQ10,'k:')
hold on
plot(Tvector(1:21),q10fMHK3,'k-') %based on MMRT-HH
plot([5 10 20],[5.3 3.3 2.6],'ok','MarkerSize',3,'MarkerEdgeColor',[0 0 0])
plot([0 10 20],[4.00 2.57 2.00],'^','MarkerSize',3,'MarkerEdgeColor',[0 0 0])
plot([0 10 20],[3.94 2.85 2.39],'sq','MarkerSize',3,'MarkerEdgeColor',[0 0 0])
plot([0 10 20],[5.25 2.88 2.32],'+','MarkerSize',3,'MarkerEdgeColor',[0 0 0])
plot([0 10 20],[7.90 3.19 2.55],'d','MarkerSize',3,'MarkerEdgeColor',[0 0 0]);
ylim([0 9])
xlabel('T ref (C)')
ylabel('Q10 fall')
box off

%% Figure 5
clear
homedir=pwd; % this should be the directory where the .m files are located
cd(homedir)
alldirs={'m01','m14','m15','m12','m03','m07','m08','m09','m10','m13','m11','m02'};
cd(homedir)
clear megaMMRT megaQ23 
c=1;
for a=1:length(alldirs)
    cd(homedir)
    [dataMMRT,dataQ23,tVec]=givemeh5Data(alldirs{a},homedir);
    megaMMRT(:,:,c)=dataMMRT;
    megaQ23(:,:,c)=dataQ23;
    c=c+1;
end

figure(20)
clf
t = tiledlayout(13,3,'TileSpacing','None');
c=1
tmax=2000
for a=1:length(alldirs)
    cd(homedir)
    
    dataMMRT=megaMMRT(:,:,c);
    dataQ23=megaQ23(:,:,c);
    
    nexttile(3*c-2)
    t2plot=(tVec>=450).*(tVec<=tmax);
    plot(tVec(logical(t2plot)),dataMMRT(1,logical(t2plot)),'k');
    hold on
    plot(tVec(2:end),dataQ23(1,:),'Color',[0.8 0.8 0.8])
    xlim([500 tmax])
    ylim([-110 50])
    set(gca,'Clipping','on')
    box off
    h = gca;
    if a ~= length(alldirs), h.XAxis.Visible = 'off'; end
    %axis off
    %axis square
    
    nexttile(3*c-1)
    plot(tVec(2:end),dataMMRT(2,:),'k');
    hold on
    plot(tVec(2:end),dataQ23(2,:),'Color',[0.8 0.8 0.8])
    xlim([500 tmax])
    ylim([-110 50])
    box off
    h = gca;
    h.YAxis.Visible = 'off';
    if a ~= length(alldirs), h.XAxis.Visible = 'off'; end
    %axis off
    %axis square
    
    %first spike
    dt=tVec(2)-tVec(1);
    t1=5;
    t2=30;
    spt=0:dt:t1+t2;
    nexttile(3*c)
    dummyMMRT=dataMMRT(1,:); %only 21 C
    spMMRT=find(diff(dummyMMRT(1,:)>0));
    dummyQ23=dataQ23(1,:);
    spQ23=find(diff(dummyQ23(1,:)>0));
    plot(spt(1:end),dummyMMRT(spMMRT(1)-t1/dt:spMMRT(1)+t2/dt),'k')
    hold on
    plot(spt(1:end),dummyQ23(spQ23(1)-t1/dt:spQ23(1)+t2/dt),'Color',[0.8 0.8 0.8])
    box off
    ylim([-110 50])
    h = gca;
    h.YAxis.Visible = 'off';
    if a ~= length(alldirs), h.XAxis.Visible = 'off'; end
    
    drawnow
    c=c+1;
end

%% Figure 6
%Based on "Dendritic sodium spikes are required for long-term potentiation at distal synapses on hippocampal pyramidal neurons"
inK=@(T) T+273.15;
R = 8.31;      % J/(mol*K)
kb = 1.38064852e-23; %m2 kg s-2 K-1
hplank = 6.62607004e-34; %m2 kg / s
T0=inK(25);


%MMRT Q10
% Na Q10: Cp=-2860  (-3791, -1930), dH=9.011e+04  (8.263e+04, 9.759e+04)
CpNa=-2860; dHNa=9.011e+04;
% K Q10:  Cp=-4144  (-5499, -2789), dH=8.898e+04  (7.716e+04, 1.008e+05)
CpK=-4144; dHK=8.898e+04;
% % Ca Q10: Cp=-3607  (-5471, -1743), dH=9.718e+04  (8.55e+04, 1.089e+05)
% CpCa=-3607; dHCa=9.718e+04;

trefNa=24;trefK=24; %mmfun(inK(tref))==1,aafun(inK(tref))==1
dSNa=58.5; dSK=54.7;

mmfunNa=@(T) (kb.*T/hplank).*exp(-(dHNa+CpNa.*(T-T0))./(R.*T)+(dSNa+CpNa.*(log(T)-log(T0)))./R);
mmfunK=@(T) (kb.*T/hplank).*exp(-(dHK+CpK.*(T-T0))./(R.*T)+(dSK+CpK.*(log(T)-log(T0)))./R);

Na35=mmfunNa(inK(35));
K35=mmfunK(inK(35));
Na21=mmfunNa(inK(21));
K21=mmfunK(inK(21));

homedir=pwd; %where the .m files are located

filename35=[homedir '/dspikeNaLTPHpc_eLife/extracteddata/VoltApicalTrunk35'];
filename21=[homedir '/dspikeNaLTPHpc_eLife/extracteddata/VoltApicalTrunk21'];
filename={filename21; filename35};

clf
sp=[121 122];
figure(6)
clf

for name=[1 2]
    subplot(sp(name))
    hold on
%     for sheet = [1 7]%[1 2 7 8]%  %  %1:6
        A(:,:,1)=xlsread(filename{name},1); %sheet=1 orig paper
        plot(A(:,1,1),A(:,2,1),'Color',[0.8 0.8 0.8])
        B(:,:,7)=xlsread(filename{name},7); %sheet=7 MMRT
        plot(B(:,1,7),B(:,2,7),'k')
%     end
    xlim([18 80])
end

subplot(121)
title('35 C')
xlabel('t (ms)')
ylabel('Voltage at apical trunk (mV)')
ylim([-80 10])
subplot(122)
title('21 C')
xlabel('t (ms)')
legend('original Q10','MMRT Q10');
legend boxoff

%% Figure 7

R = 8.31;      % J/(mol*K)
kb = 1.38064852e-23; %m2 kg s-2 K-1
hplank = 6.62607004e-34; %m2 kg / s

inK=@(T)(T+273.15);
T0=inK(25);
T2=100;
T=inK(1:T2);

%MMRT Q10
CpNa=-2860; dHNa=9.011e+04;
CpK=-4144; dHK=8.898e+04;
CpCa=-3607; dHCa=9.718e+04;

tref=20;
dSNa=63.1; dSK=59.5; dSCa=87.4;

mmfunNa=@(T) (kb.*T/hplank).*exp(-(dHNa+CpNa.*(T-T0))./(R.*T)+(dSNa+CpNa.*(log(T)-log(T0)))./R);
mmfunK=@(T) (kb.*T/hplank).*exp(-(dHK+CpK.*(T-T0))./(R.*T)+(dSK+CpK.*(log(T)-log(T0)))./R);
mmfunCa=@(T) (kb.*T/hplank).*exp(-(dHCa+CpCa.*(T-T0))./(R.*T)+(dSCa+CpCa.*(log(T)-log(T0)))./R);

mfunNa=@(T) ((T+10)./T).*exp((10.*T.*(dHNa+CpNa.*(T-T0))+50.*T.*CpNa-500.*CpNa)./(R.*T.^2.*(T+10)));
mfunK=@(T) ((T+10)./T).*exp((10.*T.*(dHK+CpK.*(T-T0))+50.*T.*CpK-500.*CpK)./(R.*T.^2.*(T+10)));
mfunCa=@(T) ((T+10)./T).*exp((10.*T.*(dHCa+CpCa.*(T-T0))+50.*T.*CpCa-500.*CpCa)./(R.*T.^2.*(T+10)));

for s=337:339
    subplot(s)
    hold on
    box off
    xlabel('T (C)')
    ylabel('rate coefficient')
    ylim([0 15])
end
subplot 337; plot(T-273.15,mmfunNa(T),'k');
subplot 338; plot(T-273.15,mmfunK(T),'k');
subplot 339; plot(T-273.15,mmfunCa(T),'k');

clear piNa piK piCa 
for dT=[1 3 10] %1:10 
    for Tfinal=20:T2
        if Tfinal>tref, trange=[tref:dT:Tfinal];
        else, trange=[tref:-dT:Tfinal];
        end
        piNa(Tfinal)=prod(mfunNa(inK(trange(1:end-1))).^(dT/10));
        piK(Tfinal)=prod(mfunK(inK(trange(1:end-1))).^(dT/10));
        piCa(Tfinal)=prod(mfunCa(inK(trange(1:end-1))).^(dT/10));
    end
    subplot 337; plot(piNa,'k--')
    subplot 338; plot(piK,'k--')
    subplot 339; plot(piCa,'k--')
end

% % % 
tbase=20;

myq10Na=mfunNa(inK(tbase)).^((0:20)./10)*mmfunNa(inK(tbase));
myq10K=mfunK(inK(tbase)).^((0:20)./10)*mmfunK(inK(tbase));
myq10Ca=mfunCa(inK(tbase)).^((0:20)./10)*mmfunCa(inK(tbase));

for s=331:333
    subplot(s)
    hold on
    box off
    xlabel('T (C)')
    ylabel('rate coefficient')
    ylim([0 15])
end
subplot 331;title('Na')
plot(20:40,mmfunNa(inK(20:40)),'k')
plot(20:40,myq10Na,'color',[.5 .5 .5])
subplot 332;title('K')
plot(20:40,mmfunK(inK(20:40)),'k')
plot(20:40,myq10K,'color',[.5 .5 .5])
subplot 333;title('Ca')
plot(20:40,mmfunCa(inK(20:40)),'k')
plot(20:40,myq10Ca,'color',[.5 .5 .5])

subplot 334
plot(20:40,100*(myq10Na-mmfunNa(inK(20:40)))./mmfunNa(inK(20:40)),'k')
subplot 335
plot(20:40,100*(myq10K-mmfunK(inK(20:40)))./mmfunK(inK(20:40)),'k')
subplot 336
plot(20:40,100*(myq10Ca-mmfunCa(inK(20:40)))./mmfunCa(inK(20:40)),'k')
for s=334:336
    subplot(s)
    hold on
    box off
    plot(20:40,zeros(1,21),'k--')
    xlabel('T (C)')
    ylabel('percent error')
    ylim([-20 100])
end
