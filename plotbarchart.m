function [fighandle] = plotbarchart(totalerr,title1,title2,vecy1,vecy2)

labelx = {'no assumption on zero location',...
    'know ~ 20% of zero location',...
    'know ~ 50% of zero location'};

modelc = {'BIC','AIC','AICC','KIC','KICC'};
vecy1str = cellstr(num2str(vecy1', '%1.2f '));
vecy2str = cellstr(num2str(vecy2', '%1.2f '));



fighandle = figure;
subplot(211); 
b = bar(totalerr(:,:,1),'FaceColor','flat'); title(title1); grid on; 
for k = 1:5
    b(k).CData = k;
end

legend(modelc,'location','southeast')
ax = gca; ax.FontSize = 20; 
if ~isempty(vecy1)
    ax.YTick = vecy1; ax.YTickLabel = vecy1str;
end
ax.XTickLabel = labelx;
set(gcf, 'Position', get(0, 'Screensize'));
[ax] = nowhitespace(ax);


subplot(212);
b = bar(totalerr(:,:,2),'FaceColor','flat'); title(title2); grid on;
for k = 1:5
    b(k).CData = k;
end

ax = gca; ax.FontSize = 20; 
if ~isempty(vecy2)
ax.YTick = vecy2; ax.YTickLabel = vecy2str;
end
ax.XTickLabel = labelx;
set(gcf, 'Position', get(0, 'Screensize'));
[ax] = nowhitespace(ax);
