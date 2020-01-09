clear
close all
%%
Ns = [1600];
transientFlags = 1;
uObs = [0, 1];
MacayealFlag = 0;
ist = [900];
T_final = [10.25, 10.5, 10.75];
dt = 0.01;
seasonType = 1;
obsT = 0.1;

%%
for i = 1: length(Ns)
    for j = 1: length(ist)
        for k = 1: length(uObs)
            for l = 1: length(T_final)
            adjointSSA(Ns(i), transientFlags(j),  uObs(k), 1- uObs(k), ...
                MacayealFlag, ist(j), T_final(l), dt, obsT, seasonType);
            end
        end
    end
end


