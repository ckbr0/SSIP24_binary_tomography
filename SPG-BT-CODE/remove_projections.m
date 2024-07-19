% Sajan, 2016.

newp=156; % new total number of projection rays;

nproj=size(M,1);


index=ones(nproj,1);

for i=1:3:nproj,
    if sum(index)>newp,
    index(i)=0;
    end;
end;

index=logical(index);

M=M(index,:);

c=M*ph_orig(:);

clear  nproj index newp i;