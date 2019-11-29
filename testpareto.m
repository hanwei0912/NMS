function pf = testpareto(name, no)

    dim = 3;
    if nargin<2, no  = 500; end
    switch lower(name)
        case {'zzj_f1','zzj_f5','zzj_f9','zzj_f10','tec09_f1','tec09_f2','tec09_f3','tec09_f4','tec09_f5','tec09_f7','tec09_f8'}
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-sqrt(pf(1,:));
        case {'zzj_f2','tec09_f9'}
            pf            = zeros(2,no);
            pf(1,:)       = linspace(0,1,no);
            pf(2,:)       = 1-pf(1,:).^2;
        case {'zzj_f3','zzj_f7'}
            pf          = zeros(2,no);
            x           = linspace(0,1,no);
            pf(1,:)     = 1-exp(-4*x).*(sin(6*pi*x)).^6;
            pf(2,:)     = 1-pf(1,:).^2;
            clear x;
        case {'zzj_f4','zzj_f8','tec09_f6'}
            weights     = pf3D(33);
            no          = size(weights,2);
            pf          = zeros(3,no);
            pf(1,:)     = cos(0.5*pi*weights(1,:)).*cos(0.5*pi*weights(2,:));
            pf(2,:)     = cos(0.5*pi*weights(1,:)).*sin(0.5*pi*weights(2,:));
            pf(3,:)     = sin(0.5*pi*weights(1,:));
            clear weights;
        case 'zzj_f6'
            pf          = zeros(2,no);
            x           = linspace(0,1,no);
            pf(1,:)     = sqrt(x);
            pf(2,:)     = 1-pf(1,:).^2;
            clear x;
        case 'uf1'
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-sqrt(pf(1,:));
        case 'uf2'
            pf            = zeros(2,no);
            pf(1,:)       = linspace(0,1,no);
            pf(2,:)       = 1-sqrt(pf(1,:));
        case 'uf3'
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-sqrt(pf(1,:));
        case 'uf4'
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-pf(1,:).^2;
        case 'uf5'
            no          = 21;
            pf          = zeros(2,no);
            pf(1,:)     = (0:1:20)/20.0;
            pf(2,:)     = 1-pf(1,:);
        case 'uf6'
            num                     = floor(no/3);
            pf                      = zeros(2,no);
            pf(1,1:num)             = 0.0;
            pf(1,(num+1):(2*num))   = linspace(0.25,0.5,num);
            pf(1,(2*num+1):no)      = linspace(0.75,1.0,no-2*num);
            pf(2,:)                 = 1-pf(1,:);
        case 'uf7'
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-pf(1,:);
        case {'uf8','uf10'}
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf          = zeros(3,no);
            pf(1,:)     = cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
            pf(2,:)     = cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));
            pf(3,:)     = sin(0.5*pi*ps(1,:));   
            clear s t;
        case 'uf9'
            num             = floor(sqrt(no));
            no              = num*num;
            noA             = floor(num/2);
            A               = zeros(1,num);
            A(1,1:noA)      = linspace(0,0.25,noA);
            A(1,noA+1:num)  = linspace(0.75,1,num-noA);
            [s,t]           = meshgrid(A,linspace(0,1,num));
            ps              = zeros(dim,no);
            ps(1,:)         = reshape(s,[1,no]);
            ps(2,:)         = reshape(t,[1,no]);            
            ps(3:dim,:)     = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf              = zeros(3,no);
            pf(1,:)         = ps(1,:).*ps(2,:);
            pf(2,:)         = (1.0-ps(1,:)).*ps(2,:);
            pf(3,:)         = 1.0-ps(2,:);    
            clear A s t;
    end
end

%%
function weights = pf3D(unit)

popsize = 0;
for i=0:1:unit, for j=0:1:unit, if i+j<=unit, popsize = popsize+1; end, end, end
weights = zeros(2, popsize);
n = 1;
for i=0:1:unit
    for j=0:1:unit
        if i+j<=unit
            weights(1,n) = i/unit;
            weights(2,n) = j/unit;
            n            = n+1;
        end
    end
end

end
