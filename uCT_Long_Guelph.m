%-----------------------------------%
% Binary image to geometry%
% This create a .geo file from a binarized image calculating boundary of a binary image and forming a surface using spline function %
% Adapted on Apr 27 2018 from (c) 2017, Sathyanarayan Rao
% @coauthor: Monique Dali
%-----------------------------------%

%File
%Sample number {'05','08'};
am = {'YZ','XZ'};
%tomograph slices
%sl = {'321','1005','1689'};
%threshold
th = ['1','2'];

for cs = am
    slice = '474';
    sc =  0.00003797; 
    for t = th
        % Insert segmented image(black-->voids, white-->porous matrix)
        image_input = strcat('Guelph_Dolomite\Pre-process\GD-',cs{1},'-',slice,'-',t,'.tif');
        text_output = strcat('gd_',cs{1},'_',slice,'_',t,'.geo');
        
        %
        [I,map] = imread(image_input);
        BW     = im2bw(I);
        BW2 = imcomplement(BW);
        Ifill = imfill(BW2,'holes');
        %Extract contour between black & white
        [B,L,N,A] = bwboundaries(BW2);
        %BW     = imbinarize(I);

        [NyI, NxI]   = size(Ifill);

        %figure; imshow(BW);
        figure;  hold on;
        % Loop through object boundaries  
        u = [];
        v = [];
        
        % Check the domains before write
%         for k = 1:N 
%             % Boundary k is the parent of a hole if the k-th column 
%             % of the adjacency matrix A contains a non-zero element 
%             if (nnz(A(:,k)) > 0) 
%                 boundary = B{k};
%                 u = [u k];
%                 c2 = boundary(:,1); 
%                 y = NyI - c2;
%                 plot(boundary(:,2),... 
%                     y,'r','LineWidth',2); 
%                 % Loop through the children of boundary k 
%                 for l = find(A(:,k))' 
%                     boundary = B{l};
%                     v = [v l];
%                     c2 = boundary(:,1); 
%                     y = NyI - c2;
%                     plot(boundary(:,2),... 
%                         y,'g','LineWidth',2); 
%                 end 
%             end 
%         end
        %saveas(gcf,graph_output)
        
        close all;
        
        cl_1         = 0.0005;
        cl_2         = 0.00025;

        % uCT scale (m)
        xi = 0.075;
        yi = 0.030;

        % domain extent
        xmin         = 1;
        xmax         = NxI; %sample size
        ymin         = 1;
        ymax         = NyI;

        %------- writing to GMsh ------------------------------------------------

        fileID = fopen(text_output,'w');

        %--- mesh size and boundary information
        fprintf(fileID,'\n');
        C  = sprintf('// mesh size description');
        fprintf(fileID,'%s\n',C);
        fprintf(fileID,'\n');
        C  = sprintf('cl_1   =  %.5f;',cl_1);
        fprintf(fileID,'%s\n',C);
        C  = sprintf('cl_2   =  %.5f;',cl_2);
        fprintf(fileID,'%s\n',C);


        fprintf(fileID,'\n');
        C  = sprintf('// boundary points');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Point(1) = {0, 0, 0, cl_1};');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Point(2) = {%.3f, 0, 0, cl_1};',xi);
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Point(3) = {%.3f, %.3f, 0, cl_1};',xi,yi);
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Point(4) = {0, %.3f, 0, cl_1};',yi);
        fprintf(fileID,'%s\n',C);
        fprintf(fileID,'\n');
        fprintf(fileID,'\n');

        C  = sprintf('// lines that connect boundary');
        fprintf(fileID,'%s\n',C);

        C  = sprintf('Line(1) = {1, 2};');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Line(2) = {3, 2};');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Line(3) = {4, 3};');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Line(4) = {1, 4};');
        fprintf(fileID,'%s\n',C);
        fprintf(fileID,'\n');

        C  = sprintf('// Mesh Parameters');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Mesh.CharacteristicLengthExtendFromBoundary = 0;');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Mesh.CharacteristicLengthMax = 0.0005;');
        fprintf(fileID,'%s\n',C);
        fprintf(fileID,'\n');

        C  = sprintf('// Define Segment coordinates');
        fprintf(fileID,'%s\n',C);

        [a,b] = size(B);

        for k = 2:a
            boundary = B{k};
            c1 = boundary(:,2);
            c2 = boundary(:,1); 
            x  = (c1-1)*sc; y = (ymax - (c2-1))*sc;
            xs = sprintf('X%d =',k);
            ys = sprintf('Y%d =',k);
            c = strjoin(arrayfun(@(x) num2str(x),x,'UniformOutput',false),',');
            d = strjoin(arrayfun(@(y) num2str(y),y,'UniformOutput',false),',');
            C = strcat(xs,'{',c,'};');
            D = strcat(ys,'{',d,'};');
            fprintf(fileID,'\n');
            fprintf(fileID,'// Hole %d\n',k);
            fprintf(fileID,'%s\n',C);
            fprintf(fileID,'%s\n',D);
        end

        fprintf(fileID,'\n');
        C  = sprintf('// Define spline surface');
        fprintf(fileID,'%s\n',C);
        fprintf(fileID,'\n');
        C  = sprintf('LN = 90;');
        fprintf(fileID,'%s\n',C);
        fprintf(fileID,'\n');

        line_loop = [];
        surfaces = [];
        i = 90;
        for k = 2:a
            j = i+1;
            m = j+1;
            C  = sprintf('nR = #X%d[ ];',k);
            fprintf(fileID,'%s\n',C);
            C  = sprintf('p0  =  newp;');
            fprintf(fileID,'%s\n',C);
            C  = sprintf('p   =  p0;');
            fprintf(fileID,'%s\n',C);
            C  = sprintf('For i In {0:nR-1}');
            fprintf(fileID,'%s\n',C);
            C  = sprintf('Point(newp)  =    {X%d[i], Y%d[i], 0, cl_2};',k,k);
            fprintf(fileID,'%s\n',C);
            C  = sprintf('EndFor');
            fprintf(fileID,'%s\n',C);
            C  = sprintf('p2  =  newp-1;');
            fprintf(fileID,'%s\n',C);
            C  = sprintf('BSpline(%i)   =  {p:p2,p};',i);
            fprintf(fileID,'%s\n',C);
            C  = sprintf('Line Loop(%i) = {%i};',j,i);
            fprintf(fileID,'%s\n',C);
            C  = sprintf('Plane Surface(%i) = {%i};',m,j);
            fprintf(fileID,'%s\n',C);
            fprintf(fileID,'\n');
            line_loop = [line_loop j];
            surfaces = [surfaces m];
            i = m+1;
        end

        l = line_loop;
        s = surfaces;

        fprintf(fileID,'\n');
        C  = sprintf('// Define all surfaces');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Line Loop(5) = {1, -2, -3, -4};');
        fprintf(fileID,'%s\n',C);
        fprintf(fileID,'\n');
        C  = sprintf('Physical Line(1) = {4};');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Physical Line(2) = {3};');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Physical Line(3) = {2};');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Physical Line(4) = {1};');
        fprintf(fileID,'%s\n',C);
        fprintf(fileID,'\n');

        xs = sprintf('Plane Surface(%i) =',m+1);
        ys = sprintf('Physical Surface(1) =');
        c = strjoin(arrayfun(@(l) num2str(l),l,'UniformOutput',false),',');
        d = strjoin(arrayfun(@(s) num2str(s),s,'UniformOutput',false),',');
        C = strcat(xs,'{5,',c,'};');
        D = strcat(ys,'{',d,'};');
        fprintf(fileID,'%s\n',C);
        C  = sprintf('Physical Surface(0) = {%i};',m+1);
        fprintf(fileID,'%s\n',C);
        fprintf(fileID,'%s\n',D);
        fclose(fileID);
        
        % open the geo file
        %uiopen(text_output,1)
    end
end
