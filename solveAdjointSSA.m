clear
close all
%%
Ns = [1600];
transientFlags = [0];
uObs = [0, 1];
MacayealFlag = 0;

%%
for i = 1: length(Ns)
    for j = 1: length(transientFlags)
        for k = 1: length(uObs)
            adjointSSA(Ns(i), transientFlags(j), uObs(k), 1- uObs(k), MacayealFlag);
        end
    end
end