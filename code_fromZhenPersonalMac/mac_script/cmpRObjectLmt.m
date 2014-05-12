function [ output_args ] = cmpRObjectLmt( numSub, numWinPerSub )
%based on window size and R object limit to estimate the maximum number of
%subjects can run the dynamic treecut
% for win69(44s): numWinPerSub=272; for win34(22s): numWinPerSub: 284
x=(numWinPerSub*4*2)^2*numSub^2-(numWinPerSub*4*2)*numSub
if x>2147483647
    output_args=0;
    disp([num2str(numSub),'subjects will exceed R object limit.'])
    else
    output_args=1;
    disp([num2str(numSub), 'subjects can be included.'])
end

end

