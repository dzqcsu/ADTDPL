function dist = adist(Xs,Xt)
    Yss = ones(size(Xs,1),1);
    Ytt = ones(size(Xt,1),1) * 2;
    
    % The results of fitclinear() may vary in a very small range, since Matlab uses SGD to optimize SVM.
    % The fluctuation is very small, ignore it
    model_linear = fitclinear([Xs;Xt],[Yss;Ytt],'learner','svm'); 
    ypred = model_linear.predict([Xs;Xt]);
    error = mae([Yss;Ytt],ypred);
    dist = 2 * (1 - 2 * error);
end