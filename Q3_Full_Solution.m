%% Setup
nx=60;                           
ny=nx/2;                                                     
dx=2/(nx-1);                     
dy=2/(ny-1);                     
x=linspace(0,2,nx);                       
y=linspace(0,1,ny);             
p=zeros(ny,nx);                  
pn=zeros(ny,nx);                 
p(:,1)=0;
p(:,nx)=y';
p(1,:)=p(2,:);                   
p(ny,:)=p(ny-1,:); 
err=1;
tol=1e-6;
k=0;
u=p;
%% Jacobi Method
while err>tol
k=k+1;
for j=2:ny-1
    for i=2:nx-1
        u(j,i)=0.25*(p(j+1,i)+p(j,i+1)+p(j-1,i)+p(j,i-1));
    end
end
u(1,:)=u(2,:);                   
u(ny,:)=u(ny-1,:);

err =sqrt(sum(sum(u-p).^2));

p=u;
p(:,1)=0;
p(:,nx)=y';
p(1,:)=u(2,:);                   
p(ny,:)=u(ny-1,:);

mid_jacobi(k)=u(ny/2,nx/2);
end
k
figure (1)
imagesc(u)
figure(2)
plot(mid_jacobi)
%% Gauss Sidel Method

nx=60;                           
ny=30;                                  
niter=10000;                     
dx=2/(nx-1);                     
dy=2/(ny-1);                     
x=linspace(0,2,nx);                       
y=linspace(0,1,ny);             
p=zeros(ny,nx);                  
pn=zeros(ny,nx);                 
p(:,1)=0;
p(:,nx)=y';
p(1,:)=p(2,:);                   
p(ny,:)=p(ny-1,:); 
err=1;
tol=1e-6;
k=0;
u=p;
while err>tol
k=k+1;
for j=2:ny-1
    for i=2:nx-1
        u(j,i)=0.25*(p(j+1,i)+u(j,i+1)+p(j-1,i)+u(j,i-1));
    end
end
u(1,:)=u(2,:);                   
u(ny,:)=u(ny-1,:);

err =sqrt(sum(sum(u-p).^2));

p=u;
p(:,1)=0;
p(:,nx)=y';
p(1,:)=u(2,:);                   
p(ny,:)=u(ny-1,:);
mid_gs(k)=u(ny/2,nx/2);
end
k
figure(3)
imagesc(u)
figure(4)
plot(mid_gs)
%% SOR


nx=60;                           
ny=30;                                                     
dx=2/(nx-1);                     
dy=2/(ny-1);                     
x=linspace(0,2,nx);                       
y=linspace(0,1,ny);             
p=zeros(ny,nx);                  
pn=zeros(ny,nx);                 
p(:,1)=0;
p(:,nx)=y';
p(1,:)=p(2,:);                   
p(ny,:)=p(ny-1,:); 
err=1;
tol=1e-6;
k=0;
u=p;
w=1.3;
while err>tol
k=k+1;
for j=2:ny-1
    for i=2:nx-1
        u(j,i)=w*0.25*(p(j+1,i)+u(j,i+1)+p(j-1,i)+u(j,i-1))+(1-w)*u(j,i);
    end
end
u(1,:)=u(2,:);                   
u(ny,:)=u(ny-1,:);

err =sqrt(sum(sum(u-p).^2));

p=u;
p(:,1)=0;
p(:,nx)=y';
p(1,:)=u(2,:);                   
p(ny,:)=u(ny-1,:);
mid_sor(k)=u(ny/2,nx/2);
end
k
figure(5)
imagesc(u)
figure(6)
plot(mid_sor)