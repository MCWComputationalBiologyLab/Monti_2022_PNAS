function varargout = Calc_vDNA(t,vDNAin0)

%% Parameters
mpar = [19063.0860579879,10.4697587913072,19.9999999985339,...
    1.33748743848498,76.3309031030219,1.77380746476150,...
    0.00100000003388195,3.10903378222424];
%vDNAin
vDNAin_kd = 0.1;

%vDNArep_max
vDNArep_max = mpar(1)*vDNAin0^2/(mpar(2)^2+vDNAin0^2);

%vDNArep_t50
vDNArep_t50 = mpar(3)*vDNAin0/(mpar(4)+vDNAin0)+mpar(5);

%nH
nH = mpar(6)*exp(-mpar(7)*vDNAin0)+mpar(8);

%% Equations
vDNAin = vDNAin0*exp(-vDNAin_kd*t);
vDNArep = vDNArep_max*t.^nH./(vDNArep_t50.^nH + t.^nH);
vDNAtot = vDNAin + vDNArep;

if (nargout == 1)
    varargout = {vDNAtot};
else
    varargout = {vDNAtot,vDNArep,vDNAin};    
end

end