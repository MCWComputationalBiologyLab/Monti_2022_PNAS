function F = FTest(E1,E2,DF1,DF2)

if (E1>E2 && DF1>DF2)
numerator = (E1-E2)/(DF1-DF2);
denominator = E2/DF2; 

F = numerator/denominator;
else %Comparing the models the wrong way or not nested models
    F = 0;
end

end