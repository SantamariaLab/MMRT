function [rise, fall,ap_amp] = apdur(spk,TIndex)
%APDUR returns the duration of AP rising and falling phase
%   The duration of the rising phase is calculated using the time points
%   where the voltage reaches an arbitrary depolarization of -45 mV and the
%   peak of the spike.
%
%   The duration of the falling phase is calculated from the peak of the
%   spike to the repolarization crossing -65 mV (in cases where the voltage
%   does not pass -65 mV during repolarization, the duration is measured
%   from the peak to 10% of the inital resting potential. 
%
%   Input: spk (time and voltage matrix -- e.g. spk(1,:) = time and
%   spk(2,:) = voltage)
%
%   Output: [rise, fall] (duration of the rising phase and falling phase of the
%   first AP given in the same units as t)
%
%
t = spk(:,1);
dt=t(2)-t(1);
v = spk(:,2)+65;
[peakV, idx] = findpeaks(v,'minpeakheight',40);
ap_amp=peakV(1);

% using dvdt for start of rise
if TIndex<30
    vrange=v(5502:end);
else
    vrange=v(5305:end);
end

if TIndex<20
    dv=diff(vrange)>0.0025+0.0025*TIndex; % /dt=*1000
else
    dv=diff(vrange)>0.05;
end

idxdV=find(dv);

if TIndex<30
    peakV10=v(idxdV(1)+5502);  
else
    peakV10=v(idxdV(1)+5305);  
end



peakV90=(peakV(1))*0.99;
% peakV10=(peakV(1))*0.20;  

vpre=v(1:idx);
riseV=(vpre>=peakV10).*(vpre<=peakV90);
riseVf=diff(find(diff([ 0; riseV; 0]))); %BP: difference between indices of Start and End of 1's
rise=riseVf.*dt;
if isempty(rise)
    rise=0;
end

fallV=(v>0).*(v<peakV90);
fallV(1:idx)=0;
diffV=find(diff([0; fallV; 0]));
fallVf=diff(diffV);
if length(fallVf)>1
    fall=fallVf(1).*dt;
    fallV(diffV(2):end)=0;
elseif length(fallVf)==1
    fall=fallVf(1).*dt;
else
    fall=0;
end
    

plot(t,v,'.k',t(logical(riseV)),v(logical(riseV)),'r',...
    t(logical(fallV)),v(logical(fallV)),'g')
ylabel('V + 65 (mV)')
xlabel('t (ms)')
title('Spike rise (50% to 90%) & fall (90% to  rest)')

end

