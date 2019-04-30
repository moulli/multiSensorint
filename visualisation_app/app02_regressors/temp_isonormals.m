PCx = rand(2000,1); PCy = rand(2000,1);
bins = 50; Msize = 10;
[N C] = hist3([PCx,PCy],[bins bins]);
CX = C{1}; CY = C{2};
N2 = N; N2(N2 == 0) = []; Nunique = unique(N2);
colors = jet(length(Nunique));
for i = 1:length(PCx)
    if isnan(PCx(i)) PCxnew(i,1) = NaN; 
        PCynew(i,1) = NaN; J(i,1) = NaN; 
    else whichoneX = find(min(abs(CX - PCx(i))) == abs(CX - PCx(i))); 
        PCxnew(i,1) = CX(whichoneX(1)); whichoneY = find(min(abs(CY - PCy(i))) == abs(CY - PCy(i))); PCynew(i,1) = CY(whichoneY(1)); J(i,1) = sub2ind([bins,bins],whichoneX(1),whichoneY(1)); 
    end
end
for i = 1:bins 
    for j = 1:bins temp = sub2ind([bins,bins],i,j); 
        Jthese = find(J == temp); 
        if ~isempty(Jthese) Ntemp = N(temp); 
            Nthis = find(Nunique == Ntemp); 
            plot(PCx(Jthese),PCy(Jthese),'.','color',colors(Nthis,:),'Markersize',Msize); 
            hold on; 
        end
    end
end
hold off;

X = [0; 1; 1; 0; 0; 1; 1; 0];
X2 = X + 1;
Y = [0; 0; 1; 1; 0; 0; 1; 1];
Z = [0; 0; 0; 0; 1; 1; 1; 1];
Xtot = [X, X2];
Ytot = [Y, Y];
Ztot = [Z, Z];
C = permute([1, 0, 0], [1, 3, 2]) .* ones(length(X), 1);
C2 = permute([0, 0, 1], [1, 3, 2]) .* ones(length(X), 1);
Ctot = [C, C2];
figure
fill3(Xtot, Ytot, Ztot, Ctot, 'FaceAlpha', 0.3)
view(3)


temp = zgrid005reg(2);
temp = gaussianize(temp, 0.01);
x = (temp.xgrid(2:end) + temp.xgrid(1:end-1)) / 2;
y = (temp.ygrid(2:end) + temp.ygrid(1:end-1)) / 2;
z = (temp.zgrid(2:end) + temp.zgrid(1:end-1)) / 2;
[X, Y, Z] = meshgrid(y, x, z);
val = zeros(size(X));
val(temp.Zindex) = temp.Zcorrel;
figure
isovalue1 = quantile(temp.Zcorrel, 0.99);
surf1 = isosurface(X, Y, Z, val, isovalue1);
p1 = patch(surf1);
isonormals(X, Y, Z, val, p1);
set(p1, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5); % set the color, mesh and transparency level of the surface
daspect([1,1,1])
view(3); axis tight
camlight; lighting gouraud
isovalue2 = quantile(temp.Zcorrel, 0.95);
surf2 = isosurface(X, Y, Z, val, isovalue2);
p2 = patch(surf2);
isonormals(X, Y, Z, val, p2);
set(p2, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
isovalue3 = quantile(temp.Zcorrel, 0.8);
surf3 = isosurface(X, Y, Z, val, isovalue3);
p3 = patch(surf3);
isonormals(X, Y, Z, val, p3);
set(p3, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.08);
isovalue4 = quantile(temp.Zcorrel, 0.1);
surf4 = isosurface(X, Y, Z, val, isovalue4);
p4 = patch(surf4);
isonormals(X, Y, Z, val, p4);
set(p4, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.04);

for i = 1:length(zgrid005reg)
    temp = char(zgrid005reg.comments(i));
    if ~isempty(regexp(temp, '1st', 'once'))
        zgrid005reg.comments(i) = string(temp(1:end-15)) + "(regressor 1)";
    elseif ~isempty(regexp(temp, '2nd', 'once'))
        zgrid005reg.comments(i) = string(temp(1:end-15)) + "(regressor 2)";
    elseif ~isempty(regexp(temp, '3rd', 'once'))
        zgrid005reg.comments(i) = string(temp(1:end-15)) + "(regressor 3)";
    elseif ~isempty(regexp(temp, '4th', 'once'))
        zgrid005reg.comments(i) = string(temp(1:end-15)) + "(regressor 4)";
    end
end






