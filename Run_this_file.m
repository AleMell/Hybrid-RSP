clc
clear all

global finalvalfx n Posvector

number_iterations = 1;

function_values = nan(number_iterations, 1);



for ii=1:number_iterations
    %initialization_noise
    initialization
    run_and_plot
    function_values(ii)=finalvalfx;
    Pos_values((ii-1)*n+1:ii*n,:)= Posvector;
    Posvector = [];
end

