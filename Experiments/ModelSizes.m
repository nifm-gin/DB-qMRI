


nn = @(p,s,z,h) (h-1)*(z.^2+z) + z*(s+1) +p*(z+1);

gllim = @(p,s,k,l) k*(1 + p*(p+1) + s*(1+p)) +s;


%%

fig = figure;

colors  = [          0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];

marks = {'o-','^-', 's-', '*-'};

P = [3 4 5 7];

Z = 100;
H = 6;
K = 50;
L = 0;


g(4) = subplot(224);
hold on

for p = 1:numel(P)
    S = 75*(p-1):300:3500;

    plot(-1,-1,marks{p}, 'Color','k')
    
    plot(S, gllim(P(p),S,K,L), marks{p}, 'Color',colors(2,:), 'HandleVisibility','off')
    plot(S, nn(P(p),S,Z,H), marks{p}, 'Color',colors(3,:), 'HandleVisibility','off');
end
ylabel('GLLiM/NN model size')  
xlabel('Signal samples (S)')
legend([repmat('P = ',numel(P),1) num2str(P')])
title('(d)')
            

S = 0:50:3000;

P = 3;
% Z = [100 200 300 500];
Z = [100 300 500 0];

marks = {'-','--', '-.', ':'};

g(1) = subplot(221);
hold on
for j = 1:length(S), s(j) = nn(P,S(j),S(j),H); end
for z = 1:numel(Z) 
    if Z(z) ~= 0 
        plot(S, nn(P,S,Z(z),H), marks{z}, 'Color',colors(3,:));
    else
        plot(S, s, marks{z}, 'Color',colors(3,:));
    end
end
ylabel('NN model size')
legend([repmat('Z = ',numel(Z),1) num2str(Z')])
title('(a)')


Z = 300;
H = [2 3 5 10];

g(2) = subplot(222);
hold on
for h = 1:numel(H)
    plot(S, nn(P,S,Z,H(h)), marks{h}, 'Color',colors(3,:));
end
legend([repmat('H = ',numel(H),1) num2str(H')])
title('(b)')



K = [50 75 100 200];

g(3) = subplot(223);
hold on
for k = 1:numel(K) 
    plot(S, gllim(P,S,K(k),L), marks{k}, 'Color',colors(2,:));
end
ylabel('GLLiM model size')  
xlabel('Signal samples (S)')
legend([repmat('K = ',numel(K),1) num2str(K')])
title('(c)')

linkaxes(g,'y')          

%%

savefig(fig, ['figures/' 'temp_model_sizes'])

