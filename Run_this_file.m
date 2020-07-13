% Runs the code for as many as "number_iterations" times. Useful when
% aiming to compare the algorithm with different random initializations.

% The vectors "function_values" and "Pos_values" store for every iteration,
% respectively, the value of the objective function at the found minimum, 
% and the whole optimization variable values at jump times.

clc
clear all

global finalvalfx n Posvector

number_iterations = 1;

function_values = nan(number_iterations, 1);


for ii=1:1:number_iterations
    initialization
    run_and_plot
    function_values(ii)=finalvalfx;
    Pos_values((ii-1)*n+1:ii*n,:)= Posvector;
    Posvector = [];
end

