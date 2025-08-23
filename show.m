close all; clear ; clc;

nx = 64; ny = nx; nz = nx;
h = 1/nx;

x = linspace(0, 4, nx);
y = linspace(0, 4, ny);
z = linspace(0, 4, nz);
[xx, yy, zz] = meshgrid(x, y, z);


figure;


    pp = 1002;  
    ss = sprintf('data1/phi%d.m', pp);
    A = load(ss);
    ss = sprintf('data2/phi%d.m', pp);
    B = load(ss);
    ss = sprintf('data3/phi%d.m', pp);
    C = load(ss);
    

    S = zeros(nx, ny, nz);
    S2 = S;
    S3 = S;
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                aa = A(ny*nz*(i-1) + nz*(j-1) + k); % nx=ny>nz
                S(i, j, k) = aa;
                aa = B(ny*nz*(i-1) + nz*(j-1) + k); % nx=ny>nz
                S2(i, j, k) = aa;
                aa = C(ny*nz*(i-1) + nz*(j-1) + k); % nx=ny>nz
                S3(i, j, k) = aa;
            end
        end
    end

    
  
    p = patch(isosurface(xx, yy, zz, S, 0));
    set(p, 'FaceColor', 'r', 'EdgeColor', 'none');hold on;
    p2 = patch(isosurface(xx, yy, zz, S2, 0));
    set(p2, 'FaceColor', [0.3010 0.7450 0.9330], 'EdgeColor', 'none');
    p3 = patch(isosurface(xx, yy, zz, S3, 0));
    set(p3, 'FaceColor', 'y', 'EdgeColor', 'none');
    daspect([1 1 1]);
    view(36, 47);
    camlight; lighting phong;
    axis([0 4 0 4 0 4]);
    box off;
    hold on;
    axis image;
    axis off;
    xlabel('x'); ylabel('y'); zlabel('z', 'rotation', 1);

tt = sprintf('cells_in_space.eps');
print('-deps', tt);
