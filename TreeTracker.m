%{
Written by Luke Weidner, Research Assistant at Colorado School of Mines
November 2020

Tree Tracker v2
Extracts, clusters, and tracks the deformation of tree trunks in point
clouds

%}

%%
t1 = readmatrix('2018-05-11.txt');
t2 = readmatrix('2018-10-12.txt');
disp('loaded clouds')
%%
threshold = 0.79; % linearity threshold

%remove points with linearity less than threshold
trunks1 = t1(t1(:,7)>threshold,1:3);
trunks2 = t2(t2(:,4)>threshold,1:3);

pc1 = pointCloud(trunks1);
pc2 = pointCloud(trunks2);

%segment both point clouds
[labs1,num1] = pcsegdist(pc1,0.5);
[labs2,num2] = pcsegdist(pc2,0.5);
disp('segmented clouds')

%sort the clusters by the number of points
count1 = histcounts(labs1,num1);
count2 = histcounts(labs2,num2);
[~,idx1] = sort(count1,'descend');
[~,idx2] = sort(count2,'descend');

%%
center1 =[];
center2 =[];
numpts1 = [];
numpts2 = [];
%for the largets N clusters (some large number), calculate the centroid
for i = 1:1000
    
    center1(i,:)=[median(trunks1(labs1==idx1(i),1)),median(trunks1(labs1==idx1(i),2)),median(trunks1(labs1==idx1(i),3))];
    center2(i,:)=[median(trunks2(labs2==idx2(i),1)),median(trunks2(labs2==idx2(i),2)),median(trunks2(labs2==idx2(i),3))];
    
end

%% Find matches and calculate deformation
numtrees = 200; % how many matches to look for
numneighbors = 3; % how many nearest neighbor clusters to look for a match

keep=NaN(numtrees,1);
count=0;
regs = {};
Results=NaN(numtrees,9);
cent1=NaN(numtrees,3);
t=zeros(numtrees,4);
centroidDifference = [];
% loop through the numtrees largest clusters in t1
for i = 1:numtrees
    iter = 1;
    match = 0;
    % find the numneighbors closest clusters in t2
    [nn,dist] = knnsearch(center2,center1(i,:),'k',numneighbors); 
    
    % start looking for a matching tree among the closest clusters. stop
    % looking once a match is found
    while iter <=numneighbors && match~=1
        
        % Create point cloud objects for the t1 and possible t2 match 
        % clusters. Use SOR to clean the clusters.
        moving = pcdenoise(pointCloud(trunks1(labs1==idx1(i),:)),'numneighbors',20);
        fixed = pcdenoise(pointCloud(trunks2(labs2==idx2(nn(iter)),:)),'numneighbors',20);
        
        % IF the two clouds have a similar number of points, AND the max
        % height is similar, AND the min height is similar, AND the
        % centroids are close, it's a match
        if (abs(length(moving.Location)-length(fixed.Location)) < 100)...
                && abs(max(moving.Location(:,3))-max(fixed.Location(:,3))) < 1.5 ...
                && abs(min(moving.Location(:,3))-min(fixed.Location(:,3))) < 1.5 ...
                && dist(iter) < 2.0

            disp(['Found a match for cloud ',num2str(i),' and the ',num2str(iter),' nearest neighbor'])
            disp(dist(iter))
            match = 1;
            
            % visualize the two clouds, manually state whether they are a
            % match or not
            figure(2)

            pcshowpair(moving,fixed)
            keep(i) = input('good? 1 = yes, 0 = no: ');
            
            if keep(i)
                % use ICP algorithm to register t1 to t2, output geometric
                % transformation matrix
                regs{i} = pcregistericp(moving,fixed,'inlierratio',0.5);
                % calculate centroid of t1 cluster (convenient for plotting
                % later)
                cent1(i,:)=median(moving.Location);
                % Use Ryan's roto-translation code to calculate the
                % displacement
                Results(i,:) = RT_Convert(moving.Location,regs{i}.T');
                count = count+1;
            end
            
            
            
            
        end
        
        iter = iter+1;
    end
    
    
    
end
%% Show deformation plot
figure(11)
XYsub = datasample(t1(:,1:2),5000);
plot(XYsub(:,1),XYsub(:,2),'.','markersize',1)
hold on
quiver(cent1(:,1),cent1(:,2),Results(:,1),Results(:,2),2)
text(cent1(:,1),cent1(:,2),num2str(Results(:,4),2))

nansum(keep)/nansum(~isnan(keep))
%% Write outputs

% writematrix([cent1,Results],'centroids.txt')


