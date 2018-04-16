classdef Pattern
    %PATTERN class with properties
    %   poly_list   : array of Polygons building the DOT pattern
    %   poly_size   : maximum expected size (circumference) of a Polygon
    %   poly_area   : maximum expected arear of a Polygon
    %   color_prop_true : true color 1 proportion w.r.t. total area
    %   color_prop_expc : expected color 1 proportion w.r.t. toal area
    %   color_balance   : default true, false if expected proportion not reached
    
    properties
        poly_list
        poly_size = 30
        poly_area = 80
        color_prop_true
        color_prop_expc
        color_balance
        cfg
    end
    
    methods
        function obj = color_polygon(obj,prop,nper,tol,dpth)
            %  assign color index to each polygon of the (global) pattern
            %
            %  with arguments:
            %   * pat  - pattern
            %   * prop - proportion       => [0:1] of polygons with color we want to generate
            %   * nper - nb of permuation => maximum permutation to find suitable color indexing
            %   * tol  - tolerance        => acceptable margin +/- tol around proportion
            %   * dpth - depth(optional)  => stop recursion when reached
            
            if nargin<5
                dpth = 10;
            end
            
            npoly = numel(obj.poly_list);
            obj.color_prop_expc = prop;
            
            for iper = 1:nper
                areas = [obj.poly_list.area];
                idx   = 1:round(prop*npoly);
                if (sum(areas(idx))/sum(areas)<(prop+tol)) && (sum(areas(idx))/sum(areas)>(prop-tol))
                    % assign colors
                    poly_colors      = 2*ones(1,npoly);
                    poly_colors(idx) = 1;
                    
                    for ipoly = 1:npoly
                        obj.poly_list(ipoly).color_index = poly_colors(ipoly);
                    end
                    trueprop = sum([obj.poly_list([obj.poly_list.color_index]==1).area])/sum([obj.poly_list.area]);
                    obj.color_prop_true = trueprop;
                    obj.color_balance = true;
                    return
                end
                obj.poly_list = obj.poly_list(randperm(npoly));
                
                if iper == nper
                    if dpth > 0
                        obj = gen_pattern;
                        obj = color_polygon(obj,prop,nper,tol,dpth-1);
                    else
                        warning('no permutation found for proportion %d\n',round(prop*100))
                        obj.color_balance = false;
                    end
                end
            end
        end
    end
    
end

