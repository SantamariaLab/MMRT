%% this funcion is used in fractional neuron integration. It integrates the fractional derivative and  the  voltage v at each time t.

function out=HHclassic(NetProp,Iinj,t,Tfinal,model2use)
load Q10z_Fidel INa IK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ncells=NetProp.Ncells;
dt=NetProp.dt;
Cm=NetProp.Cm;
v0=NetProp.v0;
vold=v0;
vrest=NetProp.vrest;
gbarK=NetProp.gK;
gbarNa=NetProp.gNa;
gbarL=NetProp.gL;
Ek=NetProp.EK;
Ena=NetProp.ENa;
El=NetProp.EL;
m=NetProp.m;
h=NetProp.h;
n=NetProp.n;
mold=m;
h1old=h;
nold=n;


Namp=NetProp.Noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vV=vrest.*ones(length(t),Ncells);




mV=m*zeros(length(t),Ncells);
hV=h*zeros(length(t),Ncells);
nV=n*zeros(length(t),Ncells);
mV(1)=m;hV(1)=h;nV(1)=n;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The weight for the  voltage memory trace of the fractional drivative for
% calculated here for the total time t for faster simulation

%Memeory weights in Caputo
%the memory trail is the sisze of memW (not infinite)



% Iionic=@(v,m,h,n,gbarL,gbarNa,gbarK,vrest,Ena,Ek)((gL*(v-EL)+...
%         gNa*m^3*h*(v-Ena)+...
%         gK*n^4*(v-Ek)));
   
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preT=1:3;%:length(t)-1;
dT=.1;
if Tfinal>6.3
    trange=[6.3:dT:Tfinal];
elseif Tfinal<6.3
    trange=[6.3:-dT:Tfinal];
else
    trange=6.3;
end

v0=-65;

inK=@(T)(T+273.15);

vdt=0.00001;
V=-100:vdt:100;

R = 8.31;      % J/(mol*K)
kb = 1.38064852e-23; %m2 kg s-2 K-1
hplank = 6.62607004e-34; %m2 kg / s

T0=inK(25);
tref=6.3;

CpNa=-2860; dHNa=9.011e+04; dSNa=84.4;
CpK=-4144; dHK=8.898e+04; dSK=83;
% CpCa=-3607; dHCa=9.718e+04;

CpNaHK=-887; dHNaHK=3.305e+04; dSNaHK=-122; %fit parameters to H&K rise
CpKHK=-2819; dHKHK=5.151e+04; dSKHK=-51; %H&K fall


mmfunNa=@(T) (kb.*T/hplank).*exp(-(dHNa+CpNa.*(T-T0))./(R.*T)+(dSNa+CpNa.*(log(T)-log(T0)))./R);
mmfunK=@(T) (kb.*T/hplank).*exp(-(dHK+CpK.*(T-T0))./(R.*T)+(dSK+CpK.*(log(T)-log(T0)))./R);

mmfunNaHK=@(T) (kb.*T/hplank).*exp(-(dHNaHK+CpNaHK.*(T-T0))./(R.*T)+(dSNaHK+CpNaHK.*(log(T)-log(T0)))./R);
mmfunKHK=@(T) (kb.*T/hplank).*exp(-(dHKHK+CpKHK.*(T-T0))./(R.*T)+(dSKHK+CpKHK.*(log(T)-log(T0)))./R);



alpham=(2.5-0.1.*(V-v0))./(exp(2.5-0.1.*(V-v0))-1);
betam=4.*exp(-(V-v0)./18);
taum=1./(alpham+betam);
ratem=1./taum;
minf=alpham./(alpham+betam);

alphah=0.07*exp(-(V-v0)/20);
betah=1./(exp(3-0.1*(V-v0))+1);
tauh=1./(alphah+betah);
rateh=1./tauh;
hinf=alphah./(alphah+betah);

alphan=(0.1-0.01*(V-v0))./(exp(1-0.1*(V-v0))-1);
betan=0.125.*exp(-(V-v0)./80);
taun=1./(alphan+betan);
raten=1./taun;
ninf=alphan./(alphan+betan);

switch upper(model2use)
    case 'MMRT-HK3'
        ratem=ratem.*mmfunNaHK(inK(Tfinal))./mmfunNaHK(inK(tref));
        rateh=rateh.*mean([mmfunKHK(inK(Tfinal))./mmfunKHK(inK(tref)),mmfunNaHK(inK(Tfinal))./mmfunNaHK(inK(tref))]);
        raten=raten.*mmfunKHK(inK(Tfinal))./mmfunKHK(inK(tref));
        
    case 'MMRT-HK2'
        ratem=ratem.*mmfunNaHK(inK(Tfinal))./mmfunNaHK(inK(tref));
        rateh=rateh.*mmfunKHK(inK(Tfinal))./mmfunKHK(inK(tref));
        raten=raten.*mmfunKHK(inK(Tfinal))./mmfunKHK(inK(tref));
        
    case 'MMRT-HK'
        ratem=ratem.*mmfunNaHK(inK(Tfinal))./mmfunNaHK(inK(tref));
        rateh=rateh.*mmfunNaHK(inK(Tfinal))./mmfunNaHK(inK(tref));
        raten=raten.*mmfunKHK(inK(Tfinal))./mmfunKHK(inK(tref));
        
    case 'MMRT-K'
        ratem=ratem.*mmfunNa(inK(Tfinal))./mmfunNa(inK(tref));
        rateh=rateh.*mmfunNa(inK(Tfinal))./mmfunNa(inK(tref));
        raten=raten.*mmfunK(inK(Tfinal))./mmfunK(inK(tref));
    case 'MMRT'
        ratem=ratem*prod(INa.m_mdl(inK(trange(1:end-1))).^(dT/10)); %dT=0 for Tfinal=6.3
        rateh=rateh*prod(INa.m_mdl(inK(trange(1:end-1))).^(dT/10));
        raten=raten*prod(IK.m_mdl(inK(trange(1:end-1))).^(dT/10));
    case 'Q10_3'
        ratem=ratem*3.^((Tfinal-6.3)/10);
        rateh=rateh*3.^((Tfinal-6.3)/10);
        raten=raten*3.^((Tfinal-6.3)/10);
    case 'ARRHENIUS'
        ratem=ratem*prod(INa.a_mdl(inK(trange(1:end-1))).^(dT/10)); %dT=0 for Tfinal=6.3
        rateh=rateh*prod(INa.a_mdl(inK(trange(1:end-1))).^(dT/10));
        raten=raten*prod(IK.a_mdl(inK(trange(1:end-1))).^(dT/10));
end
taum=1./ratem;
tauh=1./rateh;
taun=1./raten;

mf=@(Mm,p1)((minf(round(p1))-Mm)./taum(round(p1)));
hf=@(Hm,p1)((hinf(round(p1))-Hm)./tauh(round(p1)));
nf=@(Nna,p1)((ninf(round(p1))-Nna)./taun(round(p1)));

f=@(v,m,h,n,gbarL,gbarNa,gbarK,El,Ena,Ek,I,Cm)(-(1/Cm)*(gbarL*(v-El)+...
        gbarNa*m^3*h*(v-Ena)+...
        gbarK*n^4*(v-Ek))+...
        I);
for a=1:length(t)-1  
   Im=Iinj(a);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   k1=f(vold,m,h,n,gbarL,gbarNa,gbarK,El,Ena,Ek,Im,Cm);
   k2=f(vold+(dt/2)*k1,m,h,n,gbarL,gbarNa,gbarK,El,Ena,Ek,Im,Cm);
   k3=f(vold+(dt/2)*k2,m,h,n,gbarL,gbarNa,gbarK,El,Ena,Ek,Im,Cm);
   k4=f(vold+dt*k3,m,h,n,gbarL,gbarNa,gbarK,El,Ena,Ek,Im,Cm);
   v=vold+(dt/6)*(k1+2*k2+2*k3+k4);
   vV(a+1)=v;

   k1=mf(m,vold/vdt+100/vdt+1);
   k2=mf(m+(dt/2)*k1,vold/vdt+100/vdt+1);
   k3=mf(m+(dt/2)*k2,vold/vdt+100/vdt+1);
   k4=mf(m+dt*k3,vold/vdt+100/vdt+1);
   m=m+(dt/6)*(k1+2*k2+2*k3+k4);
   
   k1=hf(h,vold/vdt+100/vdt+1);
   k2=hf(h+(dt/2)*k1,vold/vdt+100/vdt+1);
   k3=hf(h+(dt/2)*k2,vold/vdt+100/vdt+1);
   k4=hf(h+dt*k3,vold/vdt+100/vdt+1);
   h=h+(dt/6)*(k1+2*k2+2*k3+k4);
   
   k1=nf(n,vold/vdt+100/vdt+1);
   k2=nf(n+(dt/2)*k1,vold/vdt+100/vdt+1);
   k3=nf(n+(dt/2)*k2,vold/vdt+100/vdt+1);
   k4=nf(n+dt*k3,vold/vdt+100/vdt+1);
   n=n+(dt/6)*(k1+2*k2+2*k3+k4);
   
    mV(a+1,:) = m;
    hV(a+1,:) = h;
    nV(a+1,:) = n;
    vold=v;

    if ~rem(a+1,round(20/dt))% if round(1/dt) >trefrac you might lose spikes
       %display(['another sec ' num2str(a*dt)])
       plot(vV)
       drawnow
    end
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

out.v=vV;
%out. VMemory= VMemory;
out.t=t;
%If you want you can save all these variables
%out.mV=mV;
%out.nV=nV;
%out.hV=hV;
%out.ratem=ratem;
%out.rateh=rateh;
%out.raten=raten;
%out.Vrange=V;
end

