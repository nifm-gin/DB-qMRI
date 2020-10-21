function [ccmap] = mycmap_extended()

depth = 25;

c1 = [0         0.1         0.4
      0         0.4470      0.7410
      0.3010    0.7450      0.9330
      0.4660    0.6740      0.1880
      0.9290    0.6940      0.1250  
      0.8500    0.3250      0.0980
      0.6350    0.0780      0.1840];

for i = 1:length(c1)-1
    g1{i} = colorGradient(c1(i,:),c1(i+1,:),depth);
end

ccmap = c1(1,:);
for i = 1:length(c1)-1
    ccmap  = [ccmap; g1{i}(2:end,:)];
end

ccmap = [ccmap; colorGradient([0.6350    0.0780      0.1840],[0 0 0],6*25)];

end