function pat = gen_pattern(cfg)
%  GEN_PATTERN  Generate DOT pattern of DOTCAT experiment
%
%  Usage: [pat] = GEN_PATTERN(cfg)
%
%  with arguments:
%   * cfg  - configuration structure cfg should contain the following fields:
%       .nang - nb of angles    => to approximate the DOT with a polygon
%       .dmtr - DOT diameter    => diameter of the DOT to draw
%       .type - type of pattern => 'size' based on polygon size, else on area
%       .dpth - depth => depth of recursion for binary space partitioning algorithm

% check input argument cfg, set defaults values
if nargin < 1
    cfg = gen_cfg;
    cfg.ntrl = nan;
    cfg.nblck = nan;
    cfg.pseudornd = nan;
    cfg.rangeStim = nan;
    cfg.rangePerc = nan;
    cfg.mean = nan;
    cfg.std = nan;
end

clear pat
%global pat?
pat = Pattern();
pat.cfg = rmfield(cfg,{'ntrl','nblck','pseudornd','rangeStim','rangePerc','mean','std'});
start_poly = Polygon();

% approximate DOT with a polygon for the binary space partioning algorithm
angl_circle = linspace(0,2*pi,cfg.nang+1);
for theta = 2:cfg.nang
    angl_circle(theta) = angl_circle(theta)+(rand-1)/80*pi;
end

for k = 1:cfg.nang
        start_poly.vertex_list(k,:) = [cfg.xy(1)/2-cos(angl_circle(k))*cfg.dmtr/2,...
                                       cfg.xy(2)/2+sin(angl_circle(k))*cfg.dmtr/2];
end

% run binary space partinioning algorithm on our approximated DOT
bsp_split(start_poly,cfg.dpth,cfg.type);

% shuffle poly_list
npoly = numel(pat.poly_list);
pat.poly_list = pat.poly_list(randperm(npoly));

    function bsp_split(polygon,depth,typ)
        % Recursive function for binary space partitioning
        %
        %  with arguments:
        %   * polygon => the area we wish to draw our pattern on
        %   * depth   => stop recursion when reached
        %   * type - type of pattern => 'size' based on polygon size, else on area
        %
        %  translated and modified from php code written by ©Ulf Åström (happyponyland.net)
        
        
        end_rec = false; % set to true to end recursion
        
        % Make a list of all edge lengths. Go through all vertex in the
        % Polygon, calculate the distance to the next to get the length of
        % the edge.
        edges = size(polygon.vertex_list,1);
        edge_len = zeros(1,edges);
        
        for i = 1:edges
            len = edge_length(polygon.vertex_list(i,:),polygon.vertex_list(mod((i),edges)+1,:));
            edge_len(i) = len;
        end
        
        % if polygon size (or polygon area) small enough, stop splitting further
        size_poly = sum(edge_len);
        area_poly = polyarea(polygon.vertex_list(:,1),polygon.vertex_list(:,2));     
        
        if strcmp(typ,'size')
            if size_poly < pat.poly_size
                end_rec = true;
            end
        else
            if area_poly < pat.poly_area
                end_rec = true;
            end
        end
        
        % If we've reached the maximum depth for our recursion.
        if depth <= 0
            end_rec = true;
        end
        
        if end_rec
            polygon.area = area_poly;
            polygon.size = size_poly;
            pat.poly_list = [pat.poly_list, polygon];
            return;
        end
        
        % This polygon needs to be divided at least once more!
        
        if strcmp(typ,'size')
            % Sort the edges by length, then pick the two longest to split
            % along. This should ensure a somewhat uniform distribution.
            [~,I] = sort(edge_len,'descend');
        else % if based on area we take 2 opposite edges
            if edges <= 4
                I = randperm(edges,2);
            else
                I = randperm(edges,2);
                while abs(diff(I))<2
                    I = randperm(edges,2);
                end
            end
%             else
%                 I(1) = randperm(edges,1);
%                 I(2) = mod(I(1)+round(edges/2)-1,edges)+1;
%             end
        end
        
        a_edge = I(1);
        b_edge = I(2);
        
        % We need to have these in order for the polygon splitting to work
        % properly. If A has a higher index than B, exchange them.
        
        if a_edge > b_edge
            % exchange
            temp = b_edge;
            b_edge = a_edge;
            a_edge = temp;
        end
        
        % These will be the new polygons.
        a_poly = Polygon();
        b_poly = Polygon();
        
        % Modulo $edges means we go back to the first point in the polygon
        % if we run over the last one.
        a1 = polygon.vertex_list(a_edge,:);
        a2 = polygon.vertex_list(mod(a_edge,edges)+1,:);
        b1 = polygon.vertex_list(b_edge,:);
        b2 = polygon.vertex_list(mod(b_edge,edges)+1,:);
        
        % Generate the dividing edge(s) C, extending from somewhere along A
        % to somewhere along B.
        c = new_edge(a1, a2, b1, b2);
        
        % Stitch together two new polygons, each consisting of "half" of the
        % old polygon. The dividing edge(s) C are included in both.
        %
        % First polygon: a_poly
        a_vertex = [polygon.vertex_list(1:a_edge,:);...
                    c;...
                    polygon.vertex_list((b_edge+1):edges,:)];
        a_poly.vertex_list = a_vertex;
        
        % Second polygon: b_poly
        c = flip(c);
        b_vertex = [polygon.vertex_list((a_edge+1):(b_edge),:);...
                    c];
        b_poly.vertex_list = b_vertex;        
        
        % Now that we have two new polygons, try to repeat the process on them.
        bsp_split(a_poly,depth-1,typ);
        bsp_split(b_poly,depth-1,typ);
    end
end

function [ret] = new_edge(a1,a2,b1,b2) %(pat,a1,a2,b1,b2)
% a1 and a2 are two Vertex of the edge A
% b1 and b2 are two Vertex of the edge B
% this function returns a new edge (defined by 2 Vertex)

a_frac = 0.4 + .2*rand;
b_frac = 0.4 + .2*rand;

% Pick a random point along the edges A and B. We will split the
% polygon between these two.
a_vert = edge_split(a1, a2, a_frac);
b_vert = edge_split(b1, b2, b_frac);

ret = [a_vert;b_vert];

end

function [V] = edge_split(a,b,fract)
%   A and B are Polygon vertex (x,y array). Returns a third point between A and
%   B. FRAC should be a number between 0 and 1 and determines the
%   distance between A - new Vertex - B (0.5 places it in the middle).

    if a(1) < b(1)
        new_x = a(1) + abs(b(1)-a(1)) * fract;
    else
        new_x = a(1) - abs(b(1)-a(1)) * fract;
    end
    if a(2) < b(2)
        new_y = a(2) + abs(b(2)-a(2)) * fract;
    else
        new_y = a(2) - abs(b(2)-a(2)) * fract;
    end
    V = [new_x,new_y];        
end

function [l] = edge_length(a,b)
%EDGE_LENGTH return the length of edge between 2 vertices a and b
    %l = hypot(b.x-a.x, b.y-a.y);
    l = hypot(b(1)-a(1), b(2)-a(2));
end
