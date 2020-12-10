function [P] = GetPorosityMap_2(V,TT,ns,hws)

    %%% Calculates porosity Map and soil volume %%%
    %%       -IN- 
    %   V:Images Stack 
    %   TT: low-upper thresholds 
    %   ns: node spacing
    %   hws: half window size
    %%      -OUT-
    %   P: Porosity map
    %   vol: total volume of sand

    intX = 1:ns:size(V,1); intY = 1:ns:size(V,2); intZ = 1:ns:size(V,3);
    P = zeros(numel(intX),numel(intY),numel(intZ),'single');

    for ix=1:numel(intX)
        for iy=1:numel(intY)
            for iz=1:numel(intZ)

                if (V(intX(ix),intY(iy),intZ(iz))==0); P(ix,iy,iz)=0 ; continue; end

                % Small cube
                C = V(max(intX(ix)-hws,1):min(intX(ix)+hws,size(V,1)),...
                        max(intY(iy)-hws,1):min(intY(iy)+hws,size(V,2)),...
                        max(intZ(iz)-hws,1):min(intZ(iz)+hws,size(V,3)));

                K = C(:) - TT; vi = sum(K>0);
                single(100*(1-(vi/sum(~(C(:)==0)))));
                P(ix,iy,iz)=single(100*(1-(vi/sum(~(C(:)==0)))));

            end
        end
    end

end