m=83;
n=86;
radius=3;
A=linspace(1,m*n,m*n);
A=reshape(A,m,n);
B=zeros(m*n,m*n);

for e=1:m*n
    q=floor(e/m);
    R=[max(mod(e,m)-radius,1):min(mod(e,m)+radius,m)];
    for i=0:radius
        idx_l=max((q-i),0)*m +R;
        idx_r=min((q+i),n)*m+R;
        B(e,idx_l)=1;
        B(e,idx_r)=1;
        
        
    end
    
    
end



% 
% for e=1:m*n
%     for i=0:r
%         l1=max(e-i*m-r,1);
%         r1=min(e-i*m+r,m*n);
%         B(e,[l1:r1])=1;
%         l2=max(e+i*m-r,1);
%         r2=min(e+i*m+r,m*n);
%         B(e,[l2:r2])=1;
%     end
% end
