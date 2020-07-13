function DeltMin = invrho(bar_ns)

global lambda_s

temp = 0;
while(rho(lambda_s*temp)<=2*bar_ns+0.0001)
   temp = temp + 0.01;
end

DeltMin = temp;
end