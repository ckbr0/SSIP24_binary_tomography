function [matrixA, matrixb]=projectionmatrix(image, angle, center, n)
% PROJECTIONMATRIX Computes the projection matrix A for discrete tomography.
%    A = PROJECTIONMATRIX(I, alpha, c, n) computes the projection matrix A for
%    the data in image I for the angle alpha centered around c. Center c
%    should be specified in image indices. n is the number of lines to the
%    left and to the right counting from the first line specified by the
%    angle alpha and center c. The total number of lines is thus 2*n+1.
%
%    [A, b] = PROJECTIONMATRIX(I, alpha, c, n) also returns the noiseless
%    projection.

% SSIP2007 - Team 16 - TP - Szeged - 2007-07-09

[nx, ny] = size(image);
num = nx * ny;

% Precompute the line directions.
a = cos(angle);
b = sin(angle);
c = - a * center(1) - b * center(2);

% Allocate matrix A and vector b.
matrixA = zeros(2*n+1, num);
matrixb = zeros(2*n+1, 1);

% We cycle through all the image pixels.
for i = 1 : nx
    for j = 1 : ny
        
        % Find the appropriate lines. Note that we can have two lines
        % passing through the same pixel if the lines are equally spaced
        % with the distance beeing the width (or height) of the voxel. The
        % approach we use is a greedy one, but the overhead is not
        % significant as we only need to check three possible lines for
        % each voxel (one is default and two are possible candidates). If
        % the distance between lines is set to be less than
        % one then we should explicitely compute all the lines beforehand.
        kt = - round(a * i + b * j + c);
        
        for k = [kt-1:1:kt+1]
            % For current line
            C = c + k;
            
            % Compute the length of the intersection. To do that we need to
            % compute the intersections between four lines that form a square
            % where the pixel (i,j) is located.
            its = [];
            L = 1;
        
            % First line. This code part is same for all four pixel
            % borders. Note that we also have a special case when the line
            % is paralel with the axes (b is zero).
            x1 = i - 0.5;
            if 0 ~= b
                y1 = (a * x1 + C) / ( -b );
                if (j - 0.5 <= y1) & (y1 <= j + 0.5)
                    its(L, :) = [x1 y1];
                    L = L + 1;
                end            
            elseif k == kt 
                its(L, :) = [x1 j];
                L = L + 1;
            end
            
            % Second line.
            x2 = i + 0.5;
            if 0 ~= b
                y2 = (a * x2 + C) / ( -b );
                if (j - 0.5 <= y2) & (y2 <= j + 0.5)
                    its(L, :) = [x2 y2];
                    L = L + 1;
                end
            elseif k == kt
                its(L, :) = [x2 j];
                L = L + 1;            
            end
            
            % Third line.
            y3 = j - 0.5;
            if 0 ~= a
                x3 = (b * y3 + C) / ( -a );
                if (i - 0.5 <= x3) & (x3 <= i + 0.5)
                    its(L, :) = [x3 y3];
                    L = L + 1;
                end
            elseif k == kt
                its(L, :) = [i y3];
                L = L + 1;
            end
            
            % Fourth line.
            y4 = j + 0.5;
            if (0 ~= a)
                x4 = (b * y4 + C) / ( -a );
                if (i - 0.5 <= x4) & (x4 <= i + 0.5)
                    its(L, :) = [x4 y4];
                    L = L + 1;
                end     
            elseif k == kt
                its(L, :) = [i y4];
                L = L + 1;            
            end
            
            % We look for the largest distance.
            d = 0;
            if (3 == L)
                d = (its(1,1) - its(2,1))^2 + (its(1,2) - its(2,2))^2;
                d = sqrt(d);
            elseif (L > 3)
                its = sortrows(its);
                d = (its(1,1) - its(L-1,1))^2 + (its(1,2) - its(L-1,2))^2;
                d = sqrt(d);
            end

            % Those are just for debbuging. Try imagesc(stored_XXX) for
            % simple visualisation.  
            stored_distances(i, j, k - kt + 2) = d;            
            
            % Store the result in the matrix A, but only if it is valid.
            % The problem could be the computation of the exact place in
            % the matrix A where we put the coefficients. That place also
            % depends on the numbering of the elements in x. MATLAB
            % stores the elements column-wise so x(:) is the vector we
            % obtain by stacking the columns together.
            if abs(k) <= abs(n)
                matrixA(k + n + 1, nx - i + 1 + (ny - j)*nx) = d;
            end
        end
        
    end     
end

% Compute the projection (quite simple once we have a matrix A).
matrixb = matrixA * image(:);