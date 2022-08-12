function ER = EvRatio(AIC1,AIC2)

delta = AIC2-AIC1;
exponent = -0.5*delta;
ER = 1/(exp(exponent));

end