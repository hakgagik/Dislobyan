function [curls, cav, div] = latticeDerivatives(centers,directions,supercell,cVect)

x=zeros(supercell);
y=zeros(supercell);
z=zeros(supercell);
u=zeros(supercell);
v=zeros(supercell);
w=zeros(supercell);

for k=1:size(centers,1)
    [k1, k2, k3]=ind2sub(supercell,k);
    x(k1,k2,k3)=centers(k,1);
    y(k1,k2,k3)=centers(k,2);
    z(k1,k2,k3)=centers(k,3);
    u(k1,k2,k3)=directions(k,1);
    v(k1,k2,k3)=directions(k,2);
    w(k1,k2,k3)=directions(k,3);
end

x(:,:,2)=x+cVect(1);
y(:,:,2)=y+cVect(2);
z(:,:,2)=z+cVect(3);
u(:,:,2)=u;
v(:,:,2)=v;
w(:,:,2)=w;

[curlX, curlY, curlZ,tempCav]=curl(x,y,z,u,v,w);
tempDiv=divergence(x,y,z,u,v,w);

curls=zeros(size(centers));
cav=zeros(size(centers,1),1);
div=zeros(size(centers,1),1);

for k1=1:supercell(1)
    for k2=1:supercell(2)
        k=sub2ind(supercell,k1,k2,1);
        curls(k,:)=[curlX(k1,k2,1),curlY(k1,k2,1),curlZ(k1,k2,1)];
        cav(k)=tempCav(k1,k2,1);
        div(k)=tempDiv(k1,k2,1);
    end
end

curlNorms=sqrt(sum(curls.^2,2));
curls((curlNorms-mean(curlNorms))>3*std(curlNorms),:)=ones(size(curls(abs(curlNorms-mean(curlNorms))>3*std(curlNorms),:)))*(mean(curlNorms)+std(curlNorms));
curlNorms=sqrt(sum(curls.^2,2));
scatter3(centers(:,1),centers(:,2),centers(:,3),450,curlNorms,'filled')
colormap cool
axis equal