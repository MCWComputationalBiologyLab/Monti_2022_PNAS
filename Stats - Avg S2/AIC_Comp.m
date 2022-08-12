function AC = AIC_Comp(AIC1,AIC2)

delta = (AIC1-AIC2);
AC = exp(-0.5*delta)/(1+exp(-0.5*delta));
%AC = 1/exp(-0.5*delta);

end