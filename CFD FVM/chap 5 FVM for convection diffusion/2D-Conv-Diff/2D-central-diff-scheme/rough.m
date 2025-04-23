% Sample plot
x = 0:0.1:10;
y = sin(x);
plot(x, y);

% Adding x-axis label
xlabel('X Axis Label', 'FontSize', 12, 'FontWeight', 'bold');  % Main x-axis label

% Create a secondary x-axis
ax = gca;  % Get current axes
ax2 = axes('Position', ax.Position, 'Color', 'none', 'XAxisLocation', 'top', 'YAxisLocation', 'right');
set(ax2, 'XLim', ax.XLim, 'YLim', [0 1]);

% Adding label to secondary x-axis
xlabel(ax2, 'x(m)', 'FontSize', 12, 'FontWeight', 'bold');  % Secondary x-axis label

% Optionally, you can also adjust the y-axis label if needed
ylabel('Y Axis Label', 'FontSize', 12, 'FontWeight', 'bold');
