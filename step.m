function step(mop,rnd)

global params population;

    % select sub-problems to evolve
    subindex = prob_select();
    nsub     = length(subindex);
    subindex = subindex(randperm(nsub));
    
    for isub = subindex
        % increase FES of each sub-problem
        params.count(isub) = params.count(isub) + params.xdim+1;
        
        % dynamic neighborhood
        if rand() < params.pns
            neighborhood = population.neighbor(1:params.xdim+1,isub);
        else
            i = randperm(params.popsize);
            neighborhood = i(1:params.xdim+1)';
        end
        
        % search by neighborhood or history 
        if rand() < rnd
            % neighborhood
            ind = de(isub, neighborhood);
            % increase FES
            params.fes  = params.fes + 1 ;
            obj     = evaluate(mop, ind);
            population.ideapoint  = min(population.ideapoint, obj);
            update(obj, ind, 1:params.popsize);
        else
            % history
            [ind,fes] = nss(isub, neighborhood,mop);
            % increase FES
            params.fes  = params.fes + fes ;
            obj     = evaluate(mop, ind);
            for k=1:fes
                population.ideapoint  = min(population.ideapoint, obj(:,k));
                update(obj(:,k), ind(:,k), 1:params.popsize);
            end
        end
    end
    
   
    % update iterate counter
    params.gen  = params.gen + 1;
    % save current population
    if params.S(1) == '2' || params.S(3) == '2' || params.S(3) == '3'
        population.hp(1:(2*params.dt-1)*params.fdim,:)       = population.hp((params.fdim+1):end,:);
        population.hp((2*params.dt-1)*params.fdim+1:end,:)   = population.objective(:,:);
    end

    % update utility
    
    clear subindex rndi ind obj neighbourhood;
end

%% update strategy
function update(obj, parameter, neighborhood)
%updating of the neighborindex with the given new individuals
global params population;

    % objective values of current solution
    newobj  = subobjective(population.W(:,neighborhood), obj, population.ideapoint, params.dmethod);
    % previous objective values
    oldobj  = subobjective(population.W(:,neighborhood), population.objective(:,neighborhood), population.ideapoint, params.dmethod);    
    % new solution is better?
    improve = (oldobj - newobj) ./ oldobj;
    % best index
    [iv,id] = max(improve);
    % a greedy replacement strategy
    if iv > 0
        population.parameter(:,id) = parameter;
        population.objective(:,id) = obj;  
        population.aindex(id)      = 1;
        params.hit(id) = params.hit(id)+1;
    end
    clear subp newobj oops oldobj C toupdate;
end

%% simplex 
function [ind,fes] = nss(index, neighborhood,mop)
global params population;
        points = population.parameter(:,neighborhood);
        points_f=population.objective(:,neighborhood);
        [~,~,ind,~,fes]     = NSS(population.W(:,index),points,1 ,2 ,0.5 ,population.ideapoint,params.name,mop.domain,points_f);
end

%%
function ind = de(index, neighborhood)

global population params;
    %parents
    si      = ones(1,3)*index;
    c       = 0;
    while si(2)==si(1) || si(3)==si(1) || si(3)==si(2)
        si(2)   = neighborhood(floor(rand()*length(neighborhood))+1);
        si(3)   = neighborhood(floor(rand()*length(neighborhood))+1);
        c       = c + 1;
        if c>10, disp('error: de'); break; end % to ensure there is no dead loop
    end

    %retrieve the individuals.
    selectpoints    = population.parameter(:, si);

    %generate new trial point
    newpoint        = selectpoints(:,1)+params.F*(selectpoints(:,2)-selectpoints(:,3));

    %repair the new value
    rnds            = rand(params.xdim,1);
    pos             = newpoint>params.xupp;
    if sum(pos)>0
        newpoint(pos) = selectpoints(pos,1) + rnds(pos,1).*(params.xupp(pos)-selectpoints(pos,1));
    end
    pos             = newpoint<params.xlow;
    if sum(pos)>0
        newpoint(pos) = selectpoints(pos,1) - rnds(pos,1).*(selectpoints(pos,1)-params.xlow(pos));
    end
    
    ind             = polymutate(newpoint, params.pm);

    clear si selectpoints newpoint pos;
end

%%
function ind = polymutate(ind, rate)
global params;

    eta_m   = params.etam;
    mut_pow = 1.0 / (eta_m + 1.0);   
    for j = 1:params.xdim
      r = rand();
      if r <= rate
        y       = ind(j);
        yl      = params.xlow(j);
        yu      = params.xupp(j);
        delta1  = (y - yl) / (yu - yl);
        delta2  = (yu - y) / (yu - yl);
        
        rnd     = rand();
        if (rnd <= 0.5) 
          xy    = 1.0 - delta1;
          val   = 2.0 * rnd + (1.0 - 2.0 * rnd) * (xy^(eta_m + 1.0));
          deltaq= (val^mut_pow) - 1.0;
        else 
          xy    = 1.0 - delta2;
          val   = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (xy^ (eta_m + 1.0));
          deltaq= 1.0 - (val^mut_pow);
        end

        y   = y + deltaq * (yu - yl);        

%         if y < yl, y = yl; end
%         if y > yu, y = yu; end
%         ind(j) = y;        
        if y >= yl && y <= yu
            ind(j) = y;
        end
      end
    end
end

%% select sub-problems to evolve
function index = prob_select()

global  params
 
if params.S(1) == '0' 
	index   = 1:params.popsize;
else
    ind     = rand(1,params.popsize) <= params.ps;
    index   = find(ind > 0);
end

% it is not likey to happen, just in case
if isempty(index)
    index   = 1:params.popsize;
end

clear ind;
end

%% utility function update
function utility_update()

global population params;

% step 1: utility function
switch params.S(1)
    case {'0','1'}  % no update
        return;
    case '2'        % improvement
        v = utility_improvement();       
    case '3'        % density in X-space
        v = utility_x_density(); 
    case '4'        % density in F-space
        v = utility_f_density();
    case '5'        % density in both X-space and F-space
        v = utility_x_density() + utility_f_density();
        v       = (v + 1.0E-50) / (max(v) + 1.0E-50);
    case '6'        % improvement + density in X-space
        v = utility_improvement() + utility_x_density();
        v       = (v + 1.0E-50) / (max(v) + 1.0E-50); 
    case '7'        % improvement + density in F-space
        v = utility_improvement() + utility_f_density();
        v       = (v + 1.0E-50) / (max(v) + 1.0E-50);
    case '8'        % improvement + density in X-F-space
        v = utility_improvement() + utility_x_density() + utility_f_density();
        v       = (v + 1.0E-50) / (max(v) + 1.0E-50);         
end

% record the utility value
params.hvs(1:(params.dt-1),:)    = params.hvs(2:end,:);
params.hvs(end,:)                = v;

% step 2: utility function aggregation
switch params.S(2)
    case '0'    % no aggreation
    case '1'    % neighborhood aggregation
        v = zeros(1, params.popsize);
        for i=1:params.popsize
            nb      = population.neighbor(1:params.niche,i);
            vt      = params.hvs(end,nb);
            v(i)    = mean(vt(:));
        end
        v = v / max(v);
        clear nb vt;        
    case '2'    % history aggregation
        v = mean(params.hvs);
        v = v / max(v);
    case '3'    % both neighborhood and history aggregation
        v = zeros(1, params.popsize);
        for i=1:params.popsize
            nb      = population.neighbor(1:params.niche,i);
            vt      = params.hvs(:,nb);
            v(i)    = mean(vt(:));
        end
        v = v / max(v);
        clear nb vt;        
end

% output
params.ps = v;
clear v;

end

%%
% improvement
function v = utility_improvement()
    global population params;

    cpop    = population.hp((2*params.dt-1)*params.fdim+1:end,:);
    hpop    = population.hp((params.dt-1)*params.fdim+1:params.dt*params.fdim,:);
    cobj    = subobjective(population.W, cpop, population.ideapoint,params.dmethod);
    hobj    = subobjective(population.W, hpop, population.ideapoint,params.dmethod);
    v       = (hobj-cobj)./hobj; 
    v       = (v + 1.0E-50) / (max(v) + 1.0E-50);
    clear cpop hpop cobj hobj;    
end
%%
% density in X-space
function v = utility_x_density()
    global population params;
    v       = zeros(1,params.popsize);
    for i=1:params.popsize
        neighbori 	= population.neighbor(1:params.niche, i);
        neighborx 	= population.parameter(:,neighbori);
        maxx        = max(neighborx,[],2); 
        minx        = min(neighborx,[],2);
        v(i)        = prod(maxx-minx);
%         v(i)        = sum(maxx-minx);
    end
    v       = (v + 1.0E-50) / (max(v) + 1.0E-50);
    clear neighbori neighborx maxx minx;        
end
%%
% density in F-space
function v = utility_f_density()
    global population params;
    v       = zeros(1,params.popsize);
    for i=1:params.popsize
        neighbori 	= population.neighbor(1:params.niche, i);
        neighborf 	= population.objective(:,neighbori);
        maxf        = max(neighborf,[],2);
        minf        = min(neighborf,[],2);
        v(i)        = prod(maxf-minf);
%         v(i)        = sum(maxx-minx);
    end
    v       = (v + 1.0E-50) / (max(v) + 1.0E-50);
    clear neighbori neighborf maxf minf;
end