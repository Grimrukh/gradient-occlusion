%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INVERSE DISTANCE WEIGHT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Vint]=IDW(xc,yc,vc,x,y,e,r1,r2,outSize,progress)
%%%%%%%%%%%%%%%%%
%%% INPUTS
%xc = stations x coordinates (columns) [vector]
%yc = stations y coordinates (rows) [vector]
%vc = variable values on the point [xc yc]
%x = interpolation points  x coordinates [vector]
%y = interpolation points y coordinates [vector]
%e = distance weight
%r1 --- 'fr' = fixed radius ;  'ng' = neighbours
%r2 --- radius lenght if r1 == 'fr' / number of neighbours if  r1 =='ng'
%outSize = square image output size
%%% OUTPUTS
%Vint --- Matrix [length(y),length(x)] with interpolated  variable values
%%% EXAMPLES
%%% --> V_spa=IDW(x1,y1,v1,x,y,-2,'ng',length(x1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simone Fatichi -- simonef@dicea.unifi.it
%   Copyright 2009
%   $Date: 2009/06/19 $ 
%   $Updated 2012/02/24 $ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 10
    progress = 0;
end
Vint=zeros(outSize);
xc=reshape(xc,1,length(xc));
yc=reshape(yc,1,length(yc));
vc=reshape(vc,1,length(vc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  strcmp(r1,'fr')
    if  (r2<=0)
        disp('Error: Radius must be positive')
        return
    end
    for i=1:length(x)
        D= sqrt((x(i)-xc).^2 +(y(i)-yc).^2);
        if min(D)==0
            Vint(y(i),x(i))=1;
            disp('Error: One or more stations have the coordinates of an interpolation point')
            return
        end
        vcc=vc(D<r2); % Eliminate all luminance values beyond radius
        D=D(D<r2); % Eliminate all distance values beyond radius
        V = vcc.*(D.^e); % Multiply luminances by distance to the power of weight
        wV = D.^e; % Get weights on their own
        if isempty(D)
            V=NaN;
        else
            V=sum(V)/sum(wV); % Sum and normalize weight sum to zero.
        end
        Vint(y(i),x(i))=V;
        if progress
            switch i
                case round(length(x)/10)
                    disp('10% shaded...');
                case round(length(x)*2/10)
                    disp('20% shaded...');
                case round(length(x)*3/10)
                    disp('30% shaded...');
                case round(length(x)*4/10)
                    disp('40% shaded...');
                case round(length(x)*5/10)
                    disp('50% shaded...');
                case round(length(x)*6/10)
                    disp('60% shaded...');
                case round(length(x)*7/10)
                    disp('70% shaded...');
                case round(length(x)*8/10)
                    disp('80% shaded...');
                case round(length(x)*9/10)
                    disp('90% shaded...');
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    if (r2 > length(vc)) || (r2<1)
        disp('Error: Number of neighbours not congruent with data')
        return
    end
    for i=1:length(x)
        for j=1:length(y)
            D= sqrt((x(i)-xc).^2 +(y(j)-yc).^2);
            if min(D)==0
                disp('Error: One or more stations have the coordinates of an interpolation point')
                return
            end
            [D,I]=sort(D);
            vcc=vc(I);
            V = vcc(1:r2).*(D(1:r2).^e);
            wV = D(1:r2).^e;
            V=sum(V)/sum(wV);
            Vint(y(j),x(i))=V;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
