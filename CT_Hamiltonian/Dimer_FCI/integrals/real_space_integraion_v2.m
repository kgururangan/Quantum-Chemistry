clear all
clc
close all

dx = 0.05; dy = 0.05; dz = 0.05;

x = -10:dx:10;
y = -10:dy:10;
z = -10:dz:10;

Nx = length(x); Ny = length(y); Nz = length(z);


XI1 = 1; XI2 = 1; RJ = 1;

p = 1; q = 2; r = 1; s = 2;

dr = 0.1;
r_grid = 0:dr:10;
Nr = length(r_grid);

leb_grid = getLebedevSphere(6);
Nang = length(leb_grid.x);

for A1 = 1:length(XI1)
    for A2 = 1:length(XI2)
        for J = 1:length(RJ)
            
            xi1 = XI1(A1); xi2 = XI2(A2);

            R1 = [0,0,0]; R2 = [RJ(J),0,0];

            Rat(1,:) = R1; Rat(2,:) = R2;
            g1.origin = R1; g2.origin = R2;
            g1.exps = xi1; g2.exps = xi2;
            g1.coeff = 1.0; g2.coeff = 1.0;
            [Smat, ~, ~, VVmat] = spatial_integrals_v2({g1,g2},Rat,[1,1]);
            fprintf('Gaussian integration gives V(%d,%d,%d,%d) = %4.6f\n',p,q,r,s,VVmat(p,q,r,s))
            
            S1 = @(x,y,z) primitive_norm(0,0,0,xi1)*...
                          exp(-xi1*( (x-R1(1))^2 + (y-R1(2))^2 + (z-R1(3))^2 ));
            S2 = @(x,y,z) primitive_norm(0,0,0,xi2)*...
                          exp(-xi2*( (x-R2(1))^2 + (y-R2(2))^2 + (z-R2(3))^2 ));
                      
            ORBS = {S1,S2};
            
            r1_grid = r_grid - norm(R1);
            r2_grid = r_grid - norm(R2);
            
            INT = 0.0;
            
            for i1 = 1:Nr
                %i1
                r1 = r1_grid(i1);
                
                w1_rad = 1.0;
                
                for j1 = 1:Nang
                    
                    x1 = r1*leb_grid.x(j1);
                    y1 = r1*leb_grid.y(j1);
                    z1 = r1*leb_grid.z(j1);
                    
                    w1_ang = leb_grid.w(j1);
                    
                    fpr = ORBS{p}(x1,y1,z1)*ORBS{r}(x1,y1,z1);
                    
                    for i2 = 1:Nr
                        r2 = r2_grid(i2);
                        
                        w2_rad = 1.0;
                        
                        if r1 ~= r2
                            V = 1/abs(r1-r2);
                        else
                            V = 0;
                        end
                        
                        for j2 = 1:Nang
                            
                            x2 = r2*leb_grid.x(j2);
                            y2 = r2*leb_grid.y(j2);
                            z2 = r2*leb_grid.z(j2);
                            
                            w2_ang = leb_grid.w(j2);
                            
                            fqs = ORBS{q}(x2,y2,z2)*ORBS{s}(x2,y2,z2);
                            
                            INT = INT + w1_rad*w2_rad*w1_ang*w2_ang*r2^2*r1^2*fpr*V*fqs;
                            %INT
                            
                        end
                    end
                end
            end
                            
             INT               
                            
                            
                            
%             % calculate real space orbital integral products
%             fqs = zeros(Nx,Ny,Nz);
%             fpr = zeros(Nx,Ny,Nz);
%             TEMP = zeros(Nx,Ny,Nz);
%             for i = 1:Nx
%                 for j = 1:Ny
%                     for k = 1:Nz
%                         fqs(i,j,k) = ORBS{q}(i,j,k)*ORBS{s}(i,j,k);
%                         fpr(i,j,k) = ORBS{p}(i,j,k)*ORBS{r}(i,j,k);
%                         TEMP(i,j,k) = fpr(i,j,k)*G(i,j,k)*fqs(i,j,k);
%                     end
%                 end
%             end
%             
%             INT1 = sum(sum(sum(TEMP,1)*dx,2)*dy,3)*dz;
%             
%             fprintf('Direct cartesian quadrature gives V(%d,%d,%d,%d) = %4.6f\n',p,q,r,s,INT1)
            
        end
    end
end

%%

clear all
clc
close all

leb_grid = getLebedevSphere(110);

[xs,ys,zs ] = sphere(50);
surf(xs,ys,zs)
hold on

scatter3(leb_grid.x, leb_grid.y, leb_grid.z, 50, 'MarkerFaceColor', [0,0,0])

