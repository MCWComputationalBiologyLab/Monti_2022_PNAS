function [Errs] = SensitivityAnalysis(mpar,Data1,Data2,Data3,vDNAin0s)
delta = 0.05;
stopPt = 0.5;
N = stopPt/delta;
Errs = cell(2,length(mpar));
GMaxes = DetermineGMaxes(mpar,vDNAin0s);
for i = 1:length(mpar)
    curr_mpar = mpar(i);
    dmpar = curr_mpar*delta;
    xPs = [];
    xMs = [];
    ErrPs = [];
    ErrMs = [];
   for j = 1:N
       mparP = mpar;
       mparP(i) = curr_mpar + j*dmpar;
       ErrP = Error_Stats(mparP,Data1,Data2,Data3,vDNAin0s,GMaxes);
       xPs = [xPs;mparP(i)];
       ErrPs = [ErrPs;ErrP];
       
       mparM = mpar;
       mparM(i) = curr_mpar - j*dmpar;
       ErrM = Error_Stats(mparM,Data1,Data2,Data3,vDNAin0s,GMaxes);
       xMs = [mparM(i);xMs];
       ErrMs = [ErrM;ErrMs];
   end
   
    x0 = mpar(i);
    Err0 = Error_Stats(mpar,Data1,Data2,Data3,vDNAin0s,GMaxes);
    xs = [xMs;x0;xPs];
    Err = [ErrMs;Err0;ErrPs];
    Errs{1,i} = xs;
    Errs{2,i} = Err;
end

end