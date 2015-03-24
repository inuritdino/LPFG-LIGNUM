function optim_df_plot(data,model)
% Simple plotting utility for the DATA and MODEL scatters fitting.
% USAGE:
%       OPTIM_DF_PLOT(DATA,MODEL)
%
% 

%% Find the number of scatters=DF's
if(isa(data,'cell'))
    N = length(data);
else
    N = 1;
end
if(isa(model,'cell'))
    if(length(model) ~= N)
        error('Error: dimensions of DATA and MODEL do not match.');
    end
else
    if(N ~= 1)
        error('Error: dimensions of DATA and MODEL do not match.');
    end
end

%%


end