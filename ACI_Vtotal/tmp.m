function tmp(ModelError)
%TMP    Create plot of datasets and fits
%   TMP(MODELERROR)
%   Creates a plot, similar to the plot in the main distribution fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1

% This function was automatically generated on 01-Oct-2013 17:20:04
 
% Data from dataset "ModelError data":
%    Y = ModelError
 
% Force all inputs to be column vectors
ModelError = ModelError(:);

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[551 179 688 475.975]);
legh_ = []; legt_ = {};   % handles and text for legend
ax_ = newplot;
set(ax_,'Box','on');
hold on;

% --- Plot data originally in dataset "ModelError data"
t_ = ~isnan(ModelError);
Data_ = ModelError(t_);
[F_,X_] = ecdf(Data_,'Function','cdf'...
              );  % compute empirical cdf
Bin_.rule = 1;
[C_,E_] = dfswitchyard('dfhistbins',Data_,[],[],Bin_,F_,X_);
[N_,C_] = ecdfhist(F_,X_,'edges',E_); % empirical pdf from cdf
h_ = bar(C_,N_,'hist');
set(h_,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
       'LineStyle','-', 'LineWidth',1);
xlabel('Data');
ylabel('Density')
legh_(end+1) = h_;
legt_{end+1} = 'ModelError data';

% Nudge axis limits beyond data limits
xlim_ = get(ax_,'XLim');
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
end

x_ = linspace(xlim_(1),xlim_(2),100);

% --- Create fit "Normal"

% Fit this distribution to get parameter values
t_ = ~isnan(ModelError);
Data_ = ModelError(t_);
% To use parameter estimates from the original fit:
%     p_ = [ 1.500341284837, 0.7701186890078];
pargs_ = cell(1,2);
[pargs_{:}] = normfit(Data_, 0.05);
p_ = [pargs_{:}];
y_ = normpdf(x_,p_(1), p_(2));
h_ = plot(x_,y_,'Color',[1 0 0],...
          'LineStyle','-', 'LineWidth',2,...
          'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_;
legt_{end+1} = 'Normal';

hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'}; 
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');
