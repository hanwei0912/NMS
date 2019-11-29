function initpub(varargin)

global params population;

    %Set up the initial setting for the MOEA/D
    %the initial population size.
    params.popsize      = 100;
    %the neighbourhood size.
    params.niche        = 10;
    %decomposition method
    params.dmethod      = 'ts';
    %the size for updating in the neighbourhood.
    params.updatesize   = 2;
    %the probability to search in neighborhood
    params.pns          = 0.9;
    %the parameter for mutation
    params.etam         = 20;
    %mutation probability
    params.pm           = 0.01;
    
    %handle the parameters, mainly about the popsize
    varargin = varargin{1};
    for k=1:2:length(varargin)
        name = varargin{k};
        value= varargin{k+1};
        switch name
            case 'problem'
                mop             = value;
            case 'popsize'
                params.popsize  = value;
            case 'niche'
                params.niche    = value;
            case 'method'
                params.dmethod  = value;
            case 'updatesize'
                params.updatesize= value;
            case 'pns'
                params.pns      = value;
            case 'etam'
                params.etam     = value;
            case 'pm'
                params.pm       = value;
        end
    end
    
    params.name     = mop.name;
    params.fdim     = mop.od;
    params.xdim     = mop.pd;
    % search domain
    params.xupp     = mop.domain(:,2);
    params.xlow     = mop.domain(:,1);
    
    % weight of subproblems
    population.W            = initweight(params.fdim, params.popsize, strcmp(params.dmethod, 'ts'));
    v                       = squareform(pdist(population.W'));
    [~, population.neighbor]= sort(v);
    params.popsize          = size(population.W,2);
    
    % initialize the population
    population.parameter	= init_pop(mop, params.popsize);
    population.objective	= evaluate(mop, population.parameter);

    % initialize the approximation model
    population.ideapoint    = min(population.objective,[],2);
    
    clear v;
end

%%
function pop = init_pop(prob, n)
%RANDOMNEW to generate n new point randomly from the mop problem given.

    if nargin==1, n = 1; end

    lowend  = prob.domain(:,1);
    span    = prob.domain(:,2)-lowend;
    pop     = rand(prob.pd, n).*(span(:,ones(1, n)))+ lowend(:,ones(1,n));

    clear lowend span point;
end