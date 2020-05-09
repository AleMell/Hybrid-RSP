function xdot = f_ex1_2_noise(x)
%--------------------------------------------------------------------------
    
global n XSt XEn TimInd LamInd Alp_jSt Alp_jEn ZInd QInd MInd PInd KInd AlphaSt AlphaEn VSt VEn DeltaInd D_jSt D_jEn Delta_jEn Delta_jSt RegDim a gamma theta taudet



% state

P = x(PInd);
Delta = x(DeltaInd);
V = x(VSt:VEn);

xdot = [a*P*Delta*V; 1; zeros(RegDim-1,1)];

end