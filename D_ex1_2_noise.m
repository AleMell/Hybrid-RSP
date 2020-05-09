function inside = D_ex1_2_noise(x) 
%--------------------------------------------------------------------------

global n XSt XEn TimInd LamInd Alp_jSt Alp_jEn ZInd QInd MInd PInd KInd AlphaSt AlphaEn VSt VEn DeltaInd D_jSt D_jEn Delta_jEn Delta_jSt RegDim a gamma theta taudet


tau=x(TimInd);

if (tau>=1)
    inside = 1;
else
    inside = 0;
end



end