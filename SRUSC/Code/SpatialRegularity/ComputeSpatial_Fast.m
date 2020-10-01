function[S] = ComputeSpatial_Fast(SpatialReg)
w=SpatialReg.Width;
h=SpatialReg.Height;
r=SpatialReg.r;
sz=[w h];
S=zeros(w*h,w*h);


for k=1:r % for each wxw matrix row
    for i=1:w % for each row in matrix w
        idx=[i:min(r,w)+i]; %
         for j=0:min(r+k-1,h-1)
            S(i + (k-1)*w , min(idx+j*w,(j+1)*w)) = 1;
        end
    
    end        
    
end

for k=r+1:h
    for i=1:w
        idx=[i:min(r,w)+i];
        for j=(k-r-1):min(r+k-1,h-1)
           S(i + (k-1)*w , min(idx+j*w,(j+1)*w)) = 1;
        end
    
    end        
    
end
S=max(S,S');

end



