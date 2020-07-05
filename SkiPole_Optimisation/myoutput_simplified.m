function stop = myoutput_simplified(x,optimvalues,state);
    stop = false;
    if isequal(state,'iter')
      results=[x,optimvalues.fval]
      dlmwrite(['Simplified_Problem/Iterative_Results/results_iteration_ps_',mat2str(optimvalues.iteration)],results)
    end
end