function [value] = C_ex1_2_noise(x) 
%--------------------------------------------------------------------------
global n XSt XEn TimInd LamInd Alp_jSt Alp_jEn ZInd QInd MInd PInd KInd AlphaSt AlphaEn VSt VEn DeltaInd D_jSt D_jEn Delta_jEn Delta_jSt RegDim a gamma theta taudet


tau = x(TimInd);

if (tau >= 0) && (tau < 1)  % flow condition
    value = 1;  % report flow
else 
    value = 0;  % do not report flow
end
end