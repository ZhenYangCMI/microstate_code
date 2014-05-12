function h=BarSpecial(data, overallWidth )
    colour = {'r','b'};
    [r,c] = size(data)
    h = zeros(c,1);
    width = overallWidth / c;
    offset = [-width/2 width/2];
    for i=1:c
        h(i) = bar(data(:,i),'FaceColor','r','BarWidth',width);   
        set(h(i),'XData',get(h(i),'XData')+offset(i));
        hold on               
    end    
end