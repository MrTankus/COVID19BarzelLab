function locs = locs_of_lattice(N)

% assume Israel is a rectangle of a x 4a
a = floor(sqrt(N/4));
b = 4*a;

areas = ones(b,a);

[y, x]= find(areas); 

locs = [x,y];

