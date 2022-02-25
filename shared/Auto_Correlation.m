function [G, r, g, dg, mask] = Auto_Correlation(I1 , mask, rmax, flag)


% contributed by Sarah Veatch.



if nargin<4, flag = 0; end  
if (nargin<3 || isempty(rmax)), rmax=100; end  
if (nargin<2 || isempty(mask)),    
    imagesc(I1); axis equal tight off;
    mask = roipoly;
end



N = sum(sum(I1.*mask)); 
A = sum(sum(mask));      

I1 = double(I1);         

L1 = size(I1, 1)+rmax; 
L2 = size(I1, 2)+rmax; 

NP = real(fftshift(ifft2(abs(fft2(mask, L1, L2)).^2))); 

%para = abs(fft2(I1.*mask,L1, L2)).^2
%Y = fft2(I1);
%imagesc(abs(fftshift(Y)))
%imagesc(abs(fftshift(para)))

G1 = A^2/N^2*real(fftshift(ifft2(abs(fft2(I1.*mask,L1, L2)).^2)))./NP; 
G = imcrop(G1, [floor(L2/2+1)-rmax, floor(L1/2+1)-rmax, 2*rmax, 2*rmax]);  


xvals = ones(1, 2*rmax+1)'*(-rmax:rmax);    
yvals = (-rmax:rmax)'*ones(1, 2*rmax+1);    
zvals = G;

[theta,r,v] = cart2pol(xvals,yvals, zvals);  

Ar = reshape(r,1, (2*rmax+1)^2);
 

Avals = reshape(v,1, (2*rmax+1)^2);

[rr,ind] = sort(Ar);                        

vv = Avals(ind);                             

r = 0:floor(max(rr));                       

[n bin] = histc(rr, r-.5);                  

for j = 1:rmax+1;                            
    m = bin==j;
    n2 = sum(m);                             
    if n2==0, vals(j)=0; er(j)=0;            
    else
        temp = m.*vv;
        %temp(isnan(temp))=0;
        g(j) = sum(temp)/n2;               
        
        temp = m.*(vv-g(j)).^2;
        %temp(isnan(temp))=0;
        dg(j) = sqrt(sum(temp))/n2; 
    end
end

r = 0:rmax;

%end

G(rmax+1, rmax+1) = 0;

if flag,
    r = 0:rmax;
    errorbar(r(2:length(r)), g(2:length(r)), dg(2:length(r)));
    axis tight
end

