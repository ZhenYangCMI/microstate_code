% This script create an atlas mask on Crad-200/179 with Prec and PCC removed accroding to Harvard-Oxford atlas

clear
clc

%Create atlas masks

atlas={'Crad179'};
N_total_ROI=179;


[Outdata1,VoxDim,Header1]=rest_readfile('/Users/zhenyang/Documents/microstate/mask/mask_used/mask_used_Resliced/Resliced_HarvardOxford-cort-maxprob-thr25-2mm_YCG.nii');
PreC_PCC=(Outdata1==300|Outdata1==301|Outdata1==310|Outdata1==311);
filename1=['/Users/zhenyang/Documents/microstate/mask/final_ROIs_',char(atlas),'/PreC_PCC.nii'];
rest_WriteNiftiImage(PreC_PCC,Header1,filename1);

[Outdata2,VoxDim,Header2]=rest_readfile(['/Users/zhenyang/Documents/microstate/mask/mask_used/mask_used_Resliced/', char(atlas),'.nii']);

OutdataReduced=Outdata2;
OutdataReduced(find(PreC_PCC))=0;
rest_WriteNiftiImage(OutdataReduced,Header2,['/Users/zhenyang/Documents/microstate/mask/final_ROIs_',char(atlas),'/reduced.nii']);

for i=1:N_total_ROI
    if length(find(OutdataReduced==i))~=length(find(Outdata2==i))
       OutdataReduced(find(Outdata2==i))=0;
    else
       OutdataReduced=OutdataReduced;
    end
end
temp=unique(OutdataReduced);
ROI_index=temp(find(temp~=0));
N_ROI=length(ROI_index); % totoal number of ROI should subtract the regions = 0
rest_WriteNiftiImage(OutdataReduced,Header2,['/Users/zhenyang/Documents/microstate/mask/final_ROIs_',char(atlas),'/final_reduced.nii']);
save(['/Users/zhenyang/Documents/microstate/mask/final_ROIs_',char(atlas),'/',char(atlas),'_ROI_index.txt'],'ROI_index', '-ascii', '-double', '-tabs')
xlswrite(['/Users/zhenyang/Documents/microstate/mask/final_ROIs_',char(atlas),'/',char(atlas),'_ROI_index.xls'],ROI_index)

disp ('Create a ROI mask on Crad-200/179 with Precuneus and PCC removed based on Havard-Oxford atlas.')

