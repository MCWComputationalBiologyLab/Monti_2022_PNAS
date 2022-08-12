function AStat = AIC(E,N,K)

AStat = N*log(E/N) + 2*K + (2*K*(K+1))/(N-K-1);

end