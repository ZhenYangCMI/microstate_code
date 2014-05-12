clear
clc
close all


dataLength='all_10min';

resultDir1=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, 'session1'];
resultDir2=['/home/data/Projects/microstate/DPARSF_preprocessed/results/645/', dataLength, filesep, 'session2'];
figDir=['/home/data/Projects/microstate/DPARSF_preprocessed/fig/645/', dataLength, '/correlTwoSessions/'];

    
    tmp1=load([resultDir1,'/avg_z_stationary_FC.mat']);
    tmp2=load([resultDir2,'/avg_z_stationary_FC.mat']);
    stationaryFC1=tmp1.avg_z_stationary_FC;
    stationaryFC1Transp=stationaryFC1';
    stationaryFC2=tmp2.avg_z_stationary_FC;
    stationaryFC2Transp=stationaryFC2';
    stationaryFCTransp=[stationaryFC1Transp,stationaryFC2Transp];
    [corStationaryFC,pValue]=corrcoef( stationaryFCTransp);
    
    figure(1)
    imagesc(corStationaryFC)
    colorbar
    title('Correlations between stationary FC')
    xlabel('Stationary FC Session 1 and Session 2')
    ylabel('Stationary FC Session 1 and Session 2')
    caxis([-1 1])
    saveas(figure(1),[figDir,'CorrelBetwTwoSessions_stationaryFC.png'])



    tmp1=load([resultDir2,'/meanZStationaryFCPartial.mat']);
    tmp2=load([resultDir2,'/avg_z_stationary_FC.mat']);
    partialFC1=tmp1.meanZStationaryFCPartial;
    partialFC1Transp=partialFC1';
    fullFC2=tmp2.avg_z_stationary_FC;
    fullFC2Transp=fullFC2';
    FCTransp=[partialFC1Transp,fullFC2Transp];
    [corStationaryFC,pValue]=corrcoef( FCTransp);
    
    figure(1)
    imagesc(corStationaryFC)
    colorbar
    title('Correlations between FC estimated with parital and full FC')
    xlabel('Partial correlation (1-4) and full correlation (5-8)')
    ylabel('Partial correlation (1-4) and full correlation (5-8)')
    caxis([-1 1])
    saveas(figure(1),[figDir,'CorrelBetwPartialandFull_sesson2.png'])