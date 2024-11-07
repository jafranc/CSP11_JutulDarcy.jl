% extrudeAndPopulate data for Jutul
% after python --
% import gsmh
% gmsh.initialize()
% gmsh.open(filename)
% gmsh.write(filename.m)

clear all
% run '/home/jfranc/codes/matlab'/mrst-2023a/startup.m
% run '/home/jfranc/codes/matlab'/mrst-2023a/startup_user.m
% mrstModule add gmsh
path_to_mesh = '/opt/misc/';
mesh_name = 'spe11b_structured.m';
[G,m] = gmshToMRST([path_to_mesh,mesh_name]);

% spe11-b case
perm = [1e-16, 1e-13, 2e-13, 5e-13, 1e-12, 2e-12, 0];
poro = [.1, .2, .2, .2, .25, .35, .0];

vperm = zeros(G.cells.num,3);
vporo = zeros(G.cells.num,1);
for i=1:G.cells.num
    vperm(i,:) = repmat(perm(G.cells.tag(i)),[3,1]);
    vporo(i) = poro(G.cells.tag(i));
end


%%
buffermult = zeros(1);
% tag fractionInA(B,C) from nodes.pos
% getWell cells


%make rock
% read poro and perm from vtk // ow extrude with python gmsh and data
rock = struct('perm' , vperm,...
         'poro', vporo,...
         'regions' , struct( 'saturation' , G.cells.tag ),...
         'bufferMult', buffermult );

%extrude G
% bufferCells and bufferFaces and BufferMult (again)
% quid de minOrgVol/maxOrgVol ?? 

G.cells.indexMap = 1:G.cells.num;
%rotate y<->z
G.nodes.coords(:,[1 2 3]) = G.nodes.coords(:,[1 3 2]);
G.cells.centroids(:,[1 2 3]) = G.cells.centroids(:,[1 3 2]);

%rotate and translate
G.nodes.coords(:,3) = - G.nodes.coords(:,3) + 1200;
G.cells.centroids(:,3) = - G.cells.centroids(:,3) + 1200;

%dispatch cases on bbox
Gmax =  max(G.nodes.coords,[],1);
Gmin =  min(G.nodes.coords,[],1);

%%
%well cells
%boxes defs
boxes = struct( 'A', [1.1,0.0,1.2 - 0.6,2.8,0.01,1.2 - 0.0],...
                'B', [0.0,0.0,1.2 - 1.2,1.1,0.01,1.2 - 0.6],...
                'C', [1.1,0.0,1.2 - 0.4,2.6,0.01,1.2 - 0.1]);
wells = struct( 'W1', [0.9,0.005, 1.2 - 0.3],...
                'W2', [1.7,0.005,1.2 - 0.7] );

name_tag = 'spe11';
if Gmax(1)>2.8
    fn = fieldnames(boxes);
    for k=1:numel(fn)
        if( isnumeric(boxes.(fn{k})) )
        % do stuff
            boxes.(fn{k}) = boxes.(fn{k}).*[3000,200,1000,3000,200,1000];
        end
    end
       fn = fieldnames(wells);
    for k=1:numel(fn)
        if( isnumeric(wells.(fn{k})) )
        % do stuff
            wells.(fn{k}) = wells.(fn{k}).*[3000,200,1000];
        end
    end

    if Gmax(3)>1
        %should be C
        name_tag = [ name_tag, 'c_', num2str(G.cells.num)];
    else
        %should be B
        name_tag = [ name_tag, 'b_', num2str(G.cells.num)];
    end
else %else should be A
    name_tag = [ name_tag, 'a_', num2str(G.cells.num)];
end
G.cells.wellCells = [find_closest(G,wells.W1), find_closest(G,wells.W2)];
G.cells.fractionInA = ( all(G.cells.centroids>boxes.A(1:3),2) .* all(G.cells.centroids<boxes.A(4:end),2));
G.cells.fractionInB = ( all(G.cells.centroids>boxes.B(1:3),2) .* all(G.cells.centroids<boxes.B(4:end),2));
G.cells.fractionInC = ( all(G.cells.centroids>boxes.C(1:3),2) .* all(G.cells.centroids<boxes.C(4:end),2));
G.cells.topCells = find(G.cells.centroids(:,3) > max(G.cells.centroids(:,3))-eps);
G.cells.botCells = find(G.cells.centroids(:,3) < min(G.cells.centroids(:,3)) + eps);

G.bufferCells = unique([ find(G.cells.centroids(:,1) == min(G.cells.centroids(:,1)));...
                        find(G.cells.centroids(:,1) == max(G.cells.centroids(:,1)))]);

G.bufferMult = ones(size(G.bufferCells)) + 5000*((G.cells.tag(G.bufferCells)>1) .* (G.cells.tag(G.bufferCells)<6));
rock.bufferMult = G.bufferMult;

% G.cartDims = [100 1 50];

% save([path_to_mesh,name_tag,'.mat'],"G","rock");
save -7 [path_to_mesh,name_tag,'.mat'] G rock


