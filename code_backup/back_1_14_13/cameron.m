subj1=zeros(4*128,156);
for i = 1:4
    for j = 1:128
        subj1((i-1)*128+j,:)=squeeze(z_Dyn_FC(i,:,j));
    end
end

cc=corrcoef(subj1');
cv=var(cc(1:128,1:128));
figure(1);imagesc(cc(1:128,1:128));colorbar();
figure(2);plot(cc(1,1:128));
figure(3);plot(cv);
