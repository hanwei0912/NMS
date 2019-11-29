function demo(problem,rnd,times)

global params maxfes

path('../problem',path); 
path('../problem/cec09',path); 
path('../public',path);

% problem = 'tec09_f1';   % change to ohter test instance

mop     = testmop(problem, 30);
popsize = 300;
maxfes  = 200000;

tic;
init('problem', mop, 'popsize', popsize, 'niche', 20, 'pns', 0.8, 'F', 0.5, 'S', '20', 'DT', 20, 'method','ts');

g = 1;
params.hit=zeros(popsize,1);

while params.fes < maxfes
    
    step(mop,rnd);
    updateplot(g);
    g =  g+1;
end
endt    = toc;

disp(endt); 

end


function savedata(name)
global population;

pareto  = population;
df      = [pareto.objective]; df = df'; 
ds      = [pareto.parameter]; ds = ds'; 
w       = [population.W];     w  = w';

save(name, 'df', 'ds','w');

clear pareto df ds w;
end
%%
function updateplot(gen)

global population params maxfes

df      = population.objective; df = df';
ds      = population.parameter; ds = ds';
len     = floor(params.popsize / 3);

subplot(1,2,1);
hold off; 
if size(df,2) == 2
    plot(df(1:len,1), df(1:len,2), 'ro', 'MarkerSize',4);hold on;
    plot(df(len+1:2*len,1), df(len+1:2*len,2), 'bo', 'MarkerSize',4);hold on;
    plot(df(2*len+1:end,1), df(2*len+1:end,2), 'mo', 'MarkerSize',4);hold on;
    xlabel('f1', 'FontSize', 6);
    ylabel('f2', 'FontSize', 6);
else
    plot3(df(:,1), df(:,2), df(:,3), 'ro', 'MarkerSize',4);
    xlabel('f1', 'FontSize', 6);
    ylabel('f2', 'FontSize', 6);
    zlabel('f3', 'FontSize', 6);
end
str     = sprintf('gen=%d, fes=%f', gen, params.fes/maxfes);
title(str, 'FontSize', 8);
box on;
drawnow;

subplot(1,2,2);
hold off; 
if size(ds,2) >= 3
    plot3(ds(1:len,1), ds(1:len,2), ds(1:len,3),'ro', 'MarkerSize',4);hold on;
    plot3(ds(len+1:2*len,1), ds(len+1:2*len,2), ds(len+1:2*len,3), 'bo', 'MarkerSize',4);hold on;
    plot3(ds(2*len+1:end,1), ds(2*len+1:end,2), ds(2*len+1:end,3),'mo', 'MarkerSize',4);hold on;
    xlabel('x1', 'FontSize', 6);
    ylabel('x2', 'FontSize', 6);
    zlabel('x3', 'FontSize', 6);    
elseif size(ds,2) >= 2
    plot(ds(:,1), ds(:,2), 'ro', 'MarkerSize',4);
    xlabel('x1', 'FontSize', 6);
    ylabel('x2', 'FontSize', 6);
end
str     = sprintf('X-space');
title(str, 'FontSize', 8);
box on;
drawnow;

clear pareto df ds;
end