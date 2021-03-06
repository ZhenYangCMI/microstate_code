clear
clc
close all


dataLength='all_10min';
session='session1';
resultDir=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, session, filesep, 'clustMean'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, '/multidimScale/'];
seedList={'seed1','seed2','seed3','seed4'};
numSeed=length(seedList);
numROI=156;

if strcmp(session, 'session1')
clustNumList=[5,7,6,12,8];
else
clustNumList=[6,9,10,12,9];
end

tmp=load([resultDir,'/clusterMean_',num2str(clustNumList(1)),'clusters_',session,'_normWin.mat']);
clustMean=tmp.finalMeanWinOfClust;

tmpSeed1=load([resultDir,'/clusterMean_',num2str(clustNumList(2)),'clusters_seed1_',session,'_normWin.mat']);
tmpSeed2=load([resultDir,'/clusterMean_',num2str(clustNumList(3)),'clusters_seed2_',session,'_normWin.mat']);
tmpSeed3=load([resultDir,'/clusterMean_',num2str(clustNumList(4)),'clusters_seed3_',session,'_normWin.mat']);
tmpSeed4=load([resultDir,'/clusterMean_',num2str(clustNumList(5)),'clusters_seed4_',session,'_normWin.mat']);
    clustMeanSeed1=tmpSeed1.finalMeanWinOfClust;
     clustMeanSeed2=tmpSeed2.finalMeanWinOfClust;
      clustMeanSeed3=tmpSeed3.finalMeanWinOfClust;
       clustMeanSeed4=tmpSeed4.finalMeanWinOfClust;
clustMeanConcate=[clustMean;clustMeanSeed1;clustMeanSeed2;clustMeanSeed3;clustMeanSeed4];
clustMeanConcateTransp=clustMeanConcate';
corClusters=corrcoef(clustMeanConcateTransp);

%% euclidean distance

D=pdist(clustMeanConcate, 'euclidean');
[Y,eigvals] = cmdscale(D);
[eigvals eigvals./max(abs(eigvals))]

% the scree plot of eigvalue as a function of the number of eigvalue to
% help decide how many dimenssions to use
figure(1)
plot(1:length(eigvals),eigvals,'bo-');
axis([1,length(eigvals),min(eigvals),max(eigvals)*1.1]);
xlabel('Eigenvalue number');
ylabel('Eigenvalue');
saveas(figure(1),[figDir,'eigvals_',session, 'allSeedsAndSeedSep_euclidean.png'])

if strcmp(session, 'session1')
labels={'\leftarrowa1','\leftarrowa2','\leftarrowa3','\leftarrowa4','\leftarrowa5','\leftarrows1c1','\leftarrows1c2','\leftarrows1c3','\leftarrows1c4','\leftarrows1c5','\leftarrows1c6','\leftarrows1c7',...
   '\leftarrows2c1','\leftarrows2c2','\leftarrows2c3','\leftarrows2c4','\leftarrows2c5','\leftarrows2c6', '\leftarrows3c1','\leftarrows3c2','\leftarrows3c3','\leftarrows3c4','\leftarrows3c5',...
   '\leftarrows3c6','\leftarrows3c7','\leftarrows3c8','\leftarrows3c9','\leftarrows3c10','\leftarrows3c11','\leftarrows3c12','\leftarrows4c1','\leftarrows4c2','\leftarrows4c3','\leftarrows4c4','\leftarrows4c5','\leftarrows4c6','\leftarrows4c7','\leftarrows4c8'};
else
    labels={'\leftarrowa1','\leftarrowa2','\leftarrowa3','\leftarrowa4','\leftarrowa5','\leftarrowa6','\leftarrows1c1','\leftarrows1c2','\leftarrows1c3','\leftarrows1c4','\leftarrows1c5','\leftarrows1c6','\leftarrows1c7','\leftarrows1c8','\leftarrows1c9',...
   '\leftarrows2c1','\leftarrows2c2','\leftarrows2c3','\leftarrows2c4','\leftarrows2c5','\leftarrows2c6', '\leftarrows2c7','\leftarrows2c8','\leftarrows2c9','\leftarrows2c10','\leftarrows3c1','\leftarrows3c2','\leftarrows3c3','\leftarrows3c4','\leftarrows3c5',...
   '\leftarrows3c6','\leftarrows3c7','\leftarrows3c8','\leftarrows3c9','\leftarrows3c10','\leftarrows3c11','\leftarrows3c12','\leftarrows4c1','\leftarrows4c2','\leftarrows4c3','\leftarrows4c4','\leftarrows4c5','\leftarrows4c6','\leftarrows4c7','\leftarrows4c8','\leftarrows4c9'};
end

% plot 2D
figure(2)

if strcmp(session, 'session1')
plot(Y(1:5,1), Y(1:5,2), 'ro', 'MarkerSize', 8, 'MarkerEdgeColor','r',...
    'MarkerFaceColor','r','LineWidth',2);
hold on
plot(Y(6:12,1), Y(6:12,2), 'b+', 'MarkerSize', 8,'LineWidth',2)
    hold on
plot(Y(13:18,1), Y(13:18,2), 'gs', 'MarkerSize', 8, 'LineWidth',2)
hold on
plot(Y(19:30,1), Y(19:30,2), 'kv', 'MarkerSize', 8,'LineWidth',2)
hold on
plot(Y(31:38,1), Y(31:38,2), 'c*', 'MarkerSize', 8,'LineWidth',2)
else
plot(Y(1:6,1), Y(1:6,2), 'ro', 'MarkerSize', 8, 'MarkerEdgeColor','r',...
    'MarkerFaceColor','r','LineWidth',2);
hold on
plot(Y(7:15,1), Y(7:15,2), 'b+', 'MarkerSize', 8,'LineWidth',2)
    hold on
plot(Y(16:25,1), Y(16:25,2), 'gs', 'MarkerSize', 8,'LineWidth',2)
hold on
plot(Y(26:37,1), Y(26:37,2), 'kv', 'MarkerSize', 8,'LineWidth',2)
hold on
plot(Y(38:46,1), Y(38:46,2), 'c*', 'MarkerSize', 8,'LineWidth',2)
end
axis(max(max(abs(Y))) * [-1.1,1.1,-1.1,1.1]); axis('square');
text(Y(:,1),Y(:,2),labels,'HorizontalAlignment','left');

% add the two grey line at x=0 and y=0
if feature('HGUsingMATLABClasses')
    hx = specgraphhelper('createConstantLineUsingMATLABClasses',...
        'LineStyle','-','Color',[.7 .7 .7],'Parent',gca);
    hx.Value = 0;
else
    hx = graph2d.constantline(0,'LineStyle','-','Color',[.7 .7 .7]);
end
changedependvar(hx,'x');
if feature('HGUsingMATLABClasses')
    hy = specgraphhelper('createConstantLineUsingMATLABClasses',...
        'LineStyle','-','Color',[.7 .7 .7],'Parent',gca);
    hy.Value = 0;
else
    hy = graph2d.constantline(0,'LineStyle','-','Color',[.7 .7 .7]);
end
changedependvar(hy,'y');

saveas(figure(2),[figDir,'Distance_2D_',session, 'allSeedsAndSeedSep_euclidean.png'])

% % plot 3D
% figure(3)
% plot3(Y(:,1), Y(:,2), Y(:,3),'bo', 'MarkerSize', 8, 'MarkerEdgeColor','k',...
%     'MarkerFaceColor','g');
% box on
% view(35, 30);
% %set(gca, 'XGrid','on')
% axis(max(max(abs(Y))) * [-1.1,1.1,-1.1,1.1,-1.1,1.1]); axis('square');
% text(Y(:,1),Y(:,2),Y(:,3),labels,'HorizontalAlignment','left');
% xlabel('x axis','Rotation', -5, 'FontSize',15);
% ylabel('y axis','Rotation',35, 'FontSize',15);
% zlabel('z axis','FontSize',15);
% saveas(figure(3),[figDir,'Distance_3D_',session, '_euclidean.png'])

%% correlation
% session can be 'session1', 'session2', or 'bothSessions'

[Y,eigvals] = cmdscale(1-corMatrix);
D=squareform(1-corMatrix);

% the scree plot of eigvalue as a function of the number of eigvalue to
% help decide how many dimenssions to use
[eigvals eigvals./max(abs(eigvals))]
figure(1)
plot(1:length(eigvals),eigvals,'bo-');
axis([1,length(eigvals),min(eigvals),max(eigvals)*1.1]);
xlabel('Eigenvalue number');
ylabel('Eigenvalue');
saveas(figure(1),[figDir,'eigvalsMultidimScal_',session, '_correl.png'])

% check whether use the first 2 columns of Y is an accurate representation
% of the original distance matrix by looking at the error in the distance
% between the 2-D congiguration and the original distancce
maxrelerr = max(abs(D - pdist(Y(:,1:2)))) / max(D)

if strcmp(session, 'session1')
labels={'\leftarrowa1','\leftarrowa2','\leftarrowa3','\leftarrowa4','\leftarrowa5','\leftarrows1c1','\leftarrows1c2','\leftarrows1c3','\leftarrows1c4','\leftarrows1c5','\leftarrows1c6','\leftarrows1c7',...
   '\leftarrows2c1','\leftarrows2c2','\leftarrows2c3','\leftarrows2c4','\leftarrows2c5','\leftarrows2c6', '\leftarrows3c1','\leftarrows3c2','\leftarrows3c3','\leftarrows3c4','\leftarrows3c5',...
   '\leftarrows3c6','\leftarrows3c7','\leftarrows3c8','\leftarrows3c9','\leftarrows3c10','\leftarrows3c11','\leftarrows3c12','\leftarrows4c1','\leftarrows4c2','\leftarrows4c3','\leftarrows4c4','\leftarrows4c5','\leftarrows4c6','\leftarrows4c7','\leftarrows4c8'};
else
    labels={'\leftarrowa1','\leftarrowa2','\leftarrowa3','\leftarrowa4','\leftarrowa5','\leftarrowa6','\leftarrows1c1','\leftarrows1c2','\leftarrows1c3','\leftarrows1c4','\leftarrows1c5','\leftarrows1c6','\leftarrows1c7','\leftarrows1c8','\leftarrows1c9',...
   '\leftarrows2c1','\leftarrows2c2','\leftarrows2c3','\leftarrows2c4','\leftarrows2c5','\leftarrows2c6', '\leftarrows2c7','\leftarrows2c8','\leftarrows2c9','\leftarrows2c10','\leftarrows3c1','\leftarrows3c2','\leftarrows3c3','\leftarrows3c4','\leftarrows3c5',...
   '\leftarrows3c6','\leftarrows3c7','\leftarrows3c8','\leftarrows3c9','\leftarrows3c10','\leftarrows3c11','\leftarrows3c12','\leftarrows4c1','\leftarrows4c2','\leftarrows4c3','\leftarrows4c4','\leftarrows4c5','\leftarrows4c6','\leftarrows4c7','\leftarrows4c8','\leftarrows4c9'};
end

% plot 2D
figure(2)

if strcmp(session, 'session1')
plot(Y(1:5,1), Y(1:5,2), 'ro', 'MarkerSize', 8, 'MarkerEdgeColor','r',...
    'MarkerFaceColor','r','LineWidth',2);
hold on
plot(Y(6:12,1), Y(6:12,2), 'b+', 'MarkerSize', 8,'LineWidth',2)
    hold on
plot(Y(13:18,1), Y(13:18,2), 'gs', 'MarkerSize', 8, 'LineWidth',2)
hold on
plot(Y(19:30,1), Y(19:30,2), 'kv', 'MarkerSize', 8,'LineWidth',2)
hold on
plot(Y(31:38,1), Y(31:38,2), 'c*', 'MarkerSize', 8,'LineWidth',2)
else
plot(Y(1:6,1), Y(1:6,2), 'ro', 'MarkerSize', 8, 'MarkerEdgeColor','r',...
    'MarkerFaceColor','r','LineWidth',2);
hold on
plot(Y(7:15,1), Y(7:15,2), 'b+', 'MarkerSize', 8,'LineWidth',2)
    hold on
plot(Y(16:25,1), Y(16:25,2), 'gs', 'MarkerSize', 8,'LineWidth',2)
hold on
plot(Y(26:37,1), Y(26:37,2), 'kv', 'MarkerSize', 8,'LineWidth',2)
hold on
plot(Y(38:46,1), Y(38:46,2), 'c*', 'MarkerSize', 8,'LineWidth',2)
end
axis(max(max(abs(Y))) * [-1.1,1.1,-1.1,1.1]); axis('square');
text(Y(:,1),Y(:,2),labels,'HorizontalAlignment','left');

% add the two grey line at x=0 and y=0
if feature('HGUsingMATLABClasses')
    hx = specgraphhelper('createConstantLineUsingMATLABClasses',...
        'LineStyle','-','Color',[.7 .7 .7],'Parent',gca);
    hx.Value = 0;
else
    hx = graph2d.constantline(0,'LineStyle','-','Color',[.7 .7 .7]);
end
changedependvar(hx,'x');
if feature('HGUsingMATLABClasses')
    hy = specgraphhelper('createConstantLineUsingMATLABClasses',...
        'LineStyle','-','Color',[.7 .7 .7],'Parent',gca);
    hy.Value = 0;
else
    hy = graph2d.constantline(0,'LineStyle','-','Color',[.7 .7 .7]);
end
changedependvar(hy,'y');

saveas(figure(2),[figDir,'Distance_',session, '_correl.png'])

% plot 3D

figure(3)
plot3(Y(:,1), Y(:,2), Y(:,3),'bo', 'MarkerSize', 8, 'MarkerEdgeColor','k',...
    'MarkerFaceColor','g');
box on
view(35, 30);
%set(gca, 'XGrid','on')
axis(max(max(abs(Y))) * [-1.1,1.1,-1.1,1.1,-1.1,1.1]); axis('square');
text(Y(:,1),Y(:,2),Y(:,3),labels,'HorizontalAlignment','left');
xlabel('x axis','Rotation', -5, 'FontSize',15);
ylabel('y axis','Rotation',35, 'FontSize',15);
zlabel('z axis','FontSize',15);
saveas(figure(3),[figDir,'Distance_3D_',session, '_correl.png'])










