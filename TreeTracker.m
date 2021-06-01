%{
Written by Luke Weidner, Research Assistant at Colorado School of Mines
November 2020

Tree Tracker v2
Extracts, clusters, and tracks the deformation of tree trunks in point
clouds

%}
tic
%% point clouds to load
t1 = readmatrix('2018-05-11new2.txt');
t2 = readmatrix('2018-10-12new2.txt');
disp('loaded clouds')
%% main parameters

Tdist=2; %max total distance between clusters
rmsethresh=0.05; %max ICP rms error for match
Zdiff=2.5; % max difference in cluster height for match
minheight = 2; % min cluster height for match
thresh1 = [0.55,0.6,0.58,170]; %thresholds: [linear 1, linear 0.5, vertical 0.5, density 0.25]
thresh2 = [0.55,0.6,0.58,200]; % the thresholds can be set differently for the two clouds [3,1,2,4] [2341]

%% advanced parameters
numneighbors = 1; % how many nearest neighbor clusters to look for a match
numtrees = 400; % how many matches to look for
nlargest = 1000; % how many largest clusters to sort to nearest neighbor pool
randomize = 0; % if randomize == 1, randomly look for numtrees matches out of nlargest clusters. if ~=1, search in order from largest cluster
LoD = 0.15; %limit of detection for plotting. Only show points with total movement larger than this value.
%% processing
label1 = TrunkThis(t1(:,4:7),thresh1,[1,2,3,4]);
label2 = TrunkThis(t2(:,4:7),thresh2,[1,2,3,4]);

%create point cloud objects
pc1 = pointCloud(t1(label1==1,1:3));
pc2 = pointCloud(t2(label2==1,1:3));
%SOR filter
pc1 = pcdenoise(pc1,'numneighbors',20);
pc2 = pcdenoise(pc2,'numneighbors',20);
%voxel downsample
pc1 = pcdownsample(pc1,'gridaverage',0.1);
pc2 = pcdownsample(pc2,'gridaverage',0.1);

%segment both point clouds
[labs1,num1] = pcsegdist(pc1,0.5);
[labs2,num2] = pcsegdist(pc2,0.5);
disp('segmented clouds')

% % plot the clustered clouds
% figure(21)
% pcshow(pc1.Location,labs1)
% colormap(prism(num1))
% set(gcf,'color','w');
% set(gca,'color','w');
% figure(22)
% pcshow(pc2.Location,labs2)
% colormap(prism(num2))

%sort the clusters by the number of points
count1 = histcounts(labs1,num1);
count2 = histcounts(labs2,num2);
[csort1,idx1] = sort(count1,'descend');
[csort2,idx2] = sort(count2,'descend');
csort1(nlargest)
csort2(nlargest)
%% Calculate centroids
center1 =[];
center2 =[];

%for the largets N clusters (some large number), calculate the centroid
for i = 1:(nlargest)
    
    center1(i,:)=[median(pc1.Location(labs1==idx1(i),1)),median(pc1.Location(labs1==idx1(i),2)),median(pc1.Location(labs1==idx1(i),3))];
    center2(i,:)=[median(pc2.Location(labs2==idx2(i),1)),median(pc2.Location(labs2==idx2(i),2)),median(pc2.Location(labs2==idx2(i),3))];
    
end

%% Find matches and calculate deformation

% keep=NaN(numtrees,1);
count=zeros(numtrees,1);
regs = {};
rmse=NaN(numtrees,1);
Results=NaN(numtrees,9);
cent1=NaN(numtrees,3);
t=zeros(numtrees,4);
centroidDifference = [];
if randomize ==1 
    randtree = randperm(nlargest,numtrees);
else
    randtree = 1:numtrees;
end
% loop through the numtrees largest clusters in t1
for i = 1:numtrees
    iter = 1;
    match = 0;
    % find the numneighbors closest clusters in t2
    [nn,dist] = knnsearch(center2,center1(randtree(i),:),'k',numneighbors);
    
    % start looking for a matching tree among the closest clusters. stop
    % looking once a match is found
    while iter <=numneighbors && match~=1
        
        % Create point cloud objects for the t1 and possible t2 match
        % clusters. Use SOR to clean the clusters.
        moving = pointCloud(pc1.Location(labs1==idx1(randtree(i)),:));
        fixed = pointCloud(pc2.Location(labs2==idx2(nn(iter)),:));
        
        
        [tempreg,~,temprmse] = pcregistericp(moving,fixed,'inlierratio',0.5);
        % Use Ryan's roto-translation code to calculate the
        % displacement
        tempresults = RT_Convert(moving.Location,tempreg.T');
        
        %code for validation and testing
%         figure(2)
%         pcshowpair(moving,fixed) 
%         set(gcf,'color','w');
%         set(gca,'color','w');
%         set(gca,'XTick',[], 'YTick', [])
%         input('ok') %to pause code to look at graph
%         tempresults(4)
%         abs(max(moving.Location(:,3))-max(fixed.Location(:,3)))
%         abs(min(moving.Location(:,3))-min(fixed.Location(:,3)))
%         abs(max(moving.Location(:,3))-min(moving.Location(:,3)))
%         abs(max(fixed.Location(:,3))-min(fixed.Location(:,3)))
%         temprmse
        
        if temprmse < rmsethresh ...
                && tempresults(4) < Tdist ...
                && abs(max(moving.Location(:,3))-max(fixed.Location(:,3))) < Zdiff ...
                && abs(min(moving.Location(:,3))-min(fixed.Location(:,3))) < Zdiff ...
                && abs(max(moving.Location(:,3))-min(moving.Location(:,3))) > minheight ...
                && abs(max(fixed.Location(:,3))-min(fixed.Location(:,3))) > minheight
            
            
            match = 1;
            disp(['Found a match for cloud ',num2str(i),' and the ',num2str(iter),' nearest neighbor'])
            
            %code for validation and testing
%             figure(2)
%             pcshowpair(moving,fixed) 
%             set(gcf,'color','w');
%             set(gca,'color','w');
%             set(gca,'XTick',[], 'YTick', [])
%             input('ok') %to pause code to look at graph
            
            
            regs{i} = tempreg;
            rmse(i) = temprmse;
            % calculate centroid of t1 cluster (convenient for plotting
            % later)
            cent1(i,:)=median(moving.Location);
            
            Results(i,:) = tempresults;
            count(i) = 1;
           
        end
        
%                     pcshowpair(pctransform(moving, regs{i}),fixed) %to observe the ICP registration results
%                     input('ok')  
        

        
        iter = iter+1;
    end
    
    
    
end
disp(['Finished. Found ',num2str(sum(count)),' matches.'])
toc
%% Show deformation plot
figure(1)
clf
sig = Results(:,4) > LoD;
pc1.Color = uint8(zeros(length(pc1.Location),3).*150);
quiver3(cent1(:,1),cent1(:,2),cent1(:,3),Results(:,1),Results(:,2),Results(:,3),2,'r')
text(cent1(:,1),cent1(:,2),cent1(:,3),num2str(Results(:,4),2),'fontsize',12)
hold on
pcshow(pc1,'markersize',1)
set(gcf,'color','w');
set(gca,'color','w');


figure(2)
dat = [cent1,Results];
direc = [sind(dat(:,9)),cosd(dat(:,9))];
scatter(dat(sig,1),dat(sig,2),40,dat(sig,7),'filled')
hold on
quiver(dat(sig,1),dat(sig,2),direc(sig,1),direc(sig,2),'k')
hold off
cb1 = colorbar;
ylabel(cb1,'Total Displacement (m)')
colormap(jet)
caxis([0,1])
xlabel('Easting (m)')
ylabel('Northing (m)')


%% Write outputs
% Attpairs{6} = [cent1,Results];
% save('Attpairs','Attpairs')
out = [cent1,Results];
out = out(~isnan(out(:,1)),:);
% writematrix(out,'centroids_attachie_19-20cropcrop2.txt')

