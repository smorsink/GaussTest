function [Fit, auxOutput] = fitness(X,extPar)
auxOutput = cell(1,size(X,2));
% This version is for comparison with Slavko's synthetic data
% This version makes use of a fixed background
for i=size(X,2):-1:1  % we want to count backwards here
    % Apply constraints
    %isPhysical=applyConstraints(X(:,i),extPar);
    disp('------ Fitness ------');
    %if isPhysical
        % Parameters from vector X.
        tic;
          
        xmean=X(1,i);
        xsigma=X(2,i);
        
        
    
        
        time=toc

        
        
        cmd = '[Fit(i),auxOutput{i}] = LogLMex(extPar.fixed.numdatabins,extPar.fixed.bins,xmean,xsigma,extPar.fixed.data)';
       
        

       %disp(cmd)
        disp(i);
        eval(cmd);
        
        time=toc
    

    
    pause(0.001);
end


