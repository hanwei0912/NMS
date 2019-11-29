function init(varargin)

global params population

    % call global init function
    initpub(varargin);
    
    % Set up the initial setting for the MOEA/D
    %the parameter for the DE.
    params.F 	= 0.5;

    %handle the parameters, mainly about the popsize        
    for k=1:2:length(varargin)
        switch varargin{k}
            case 'F' 
                params.F = varargin{k+1};
            case 'S'
                params.S = varargin{k+1};
            case 'DT'
                params.dt= varargin{k+1};     
        end
    end
    params.pm       = 1.0/params.xdim;
        
    % init utility values
    params.ps           = ones(1,params.popsize);
    params.count        = zeros(1,params.popsize);
    params.fes          = 0;
    params.gen          = 0;
   
    population.hp       = repmat(population.objective,2*params.dt,1);
    params.hvs          = zeros(params.dt, params.popsize);
    population.iextreme = sum(population.W > 0.9999) >= 1;
  
    % init utility values
    if params.S(1) == '1' || params.S(1) == '6'
        switch upper(params.name)
            case 'TEC09_F2'
                filename = 'sop/T1_6.txt';
            case 'TEC09_F3'
                filename = 'sop/T2_6.txt';
            case 'TEC09_F4'
                filename = 'sop/T3.txt';                    
        end
        p = load(filename); 
        p = p(end:-1:1);
        params.ps = p./max(p);
        if params.S(1) == '6'
            params.ps(population.iextreme) = 1.0;
        end
    end
end
