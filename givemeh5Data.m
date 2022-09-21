function [dataMMRT,dataQ23,tVec]=givemeh5Data(mdir,homedir)
basedir=[homedir '/AllenModels/'];
dir=[basedir mdir '/sims/'];
dataMMRT=[]; dataQ23=[];
for model=["MMRT/","Q23/"] 
    cd(dir)
    count=0;
    for T=[21,34,38]
        count=count+1;
        cd(dir)
        subdir=[model+"output_iclamp"+T+model];
        cd(subdir)
        if model=='MMRT/'
            dataMMRT(count,:)=h5read('v_report.h5','/report/network_name/data');
        else
            dataQ23(count,:)=h5read('v_report.h5','/report/network_name/data');
        end
    end
end

tV=h5read('v_report.h5','/report/network_name/mapping/time');
tVec=tV(1):tV(3):tV(2); %start:dt:stop
