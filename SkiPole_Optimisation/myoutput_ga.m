function stop = myoutput_nelder(x,optimvalues,state);
    stop = false;
    if isequal(state,'iter')
      results=[x,optimvalues.fval];
      dlmwrite(['Own_Optimisation/Iterative_Results/results_iteration_ga',mat2str(optimvalues.iteration)],results)
    end
end