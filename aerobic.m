clear all;
zL=130;  %Height below which [H2S] approach 0, must be large enough!
N=zL*20; %Fine grid to resolve boundary layer, 200 per meter is better, usually 2000 for saving time(also my computer)(T_T)
dz=zL/N; %grid size; 
zz=linspace(0,zL,N+1);zz=zz';

NPP = 0.1 * 1e3 / 24/60/60; % mol/m2/day -> mmol/m2/s
Rox = 1e-6 * 1e-3; % 1/mM/s, m3/mol/s -> m3/mmol/s
sink = 1e-4; %  m/s
DO2 = 1e-2 * 1e-4; % cm2/s -> m2/s
% s = 1 / 100 / 1000/365/24/60/60; % cm/kyr -> m/s

CHs = 10 .^ (-3:0.5:1);
ret = zeros(length(CHs), 2);
for CHi = 1:length(CHs)
    welling = 0;
    CH0 = CHs(CHi) * NPP / sink; % uM, mmol/m3
    oxygen0 = 150; % uM, mmol/m3

    oxygen = oxygen0 + 0 * zz;
    CH = CH0 + 0 * zz;
    
    for iter = 1:1000
        oxygenp = oxygen;
        
        Gu=[1:1:N-1 2:1:N 1:1:N-1 N];
        Gv=[2:1:N 1:1:N-1 1:1:N-1 N];
        GS=[ones(1,N-1)*(DO2/dz^2-welling/2/dz) ones(1,N-1)*(DO2/dz^2+welling/2/dz) -2*DO2/dz^2-Rox*CH(2:N)' -DO2/dz^2-welling/2/dz-Rox*CH(N+1)];
        G=sparse(Gu,Gv,GS);
        
        RHS=0*zz(1:N);
        RHS(1)=-oxygen0*(DO2/dz^2+welling/2/dz);
    
        f = RHS;
        n = length(f);
        v = zeros(n,1);   
        y = v;
        w = G(1,1);
        y(1) = f(1)/w;
        for j=2:n
            v(j-1) = G(j-1,j)/w;
            w = G(j,j) - G(j,j-1)*v(j-1);
            y(j) = ( f(j) - G(j,j-1)*y(j-1) )/w;
        end
        for j=n-1:-1:1
            y(j) = y(j) - v(j)*y(j+1);
        end
    
        oxygen = [oxygen0; y];
        
        Ioxygen = zeros(1, N+1);
        for i = 2:N+1
            Ioxygen(i) = ((oxygen(1) + oxygen(i)) / 2 + sum(oxygen(2:i-1))) * dz;
        end
        
        CHp = CH;
        for i = 1:N+1
            CH(i) = CH0 * exp(-Rox / (sink + welling) * Ioxygen(i));
        end
        if max(abs(CH - CHp)) < 1e-3 && max(abs(oxygenp - oxygen)) < 1e-3
            break;
        end
    end
    
    if find(oxygen<4.5)
        ret(CHi, 1) = (find(oxygen<4.5, 1) - 1) / 20;
        ret(CHi, 2) = CH(find(oxygen<4.5, 1));
    else
        ret(CHi, 1) = zL;
        ret(CHi, 2) = CH(end);
    end
    writematrix([zz(1:20:end) oxygen(1:20:end) CH(1:20:end)], 'Results.xlsx', 'Sheet', ['CH0=' num2str(log10(CHs(CHi))) 'NPP']);
end

writematrix([["CH0_uM", "PP_g/cm2/kyr", "penD_m", "presM_uM"]; [CHs'*NPP/sink CHs'*sink*1e-3*12/10000*60*60*24*365*1000 ret]], 'Results.xlsx', 'Sheet', 'anoxia');