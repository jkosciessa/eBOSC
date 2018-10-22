function addTaskTiming(limits)

    hold on;
    line([3 3],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-'); 
    line([3.2 3.2],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-'); 
    line([0 0],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-'); 
    line([-.2 -.2],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-')
    line([-1.2 -1.2],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-'); 
    line([-1.4 -1.4],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-')
    
end