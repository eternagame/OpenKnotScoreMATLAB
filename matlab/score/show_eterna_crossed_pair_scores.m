function show_eterna_crossed_pair_scores( plot_title, eterna_scores_twist, cc_twist, crossed_pair_scores_twist, crossed_pair_quality_scores_twist);
% show_eterna_crossed_pair_scores( plot_title, eterna_scores_twist, cc_twist, crossed_pair_scores_twist, crossed_pair_quality_scores_twist);
%
% (C) R. Das, HHMI & Stanford University, 2023.

set(figure(11),'pos',[200 200 800 800],'color','white');
subplot(2,2,1);plot(eterna_scores_twist, cc_twist,'.')
xlabel( 'Eterna classic score');
ylabel( 'SHAPE Correlation coefficient');
title( plot_title,'Interpreter','none')

subplot(2,2,2);plot(eterna_scores_twist, crossed_pair_scores_twist,'.')
xlabel( 'Eterna classic score');
ylabel( 'Crossed pair score');

subplot(2,2,3);plot(eterna_scores_twist, crossed_pair_quality_scores_twist,'.')
xlabel( 'Eterna classic score');
ylabel( 'Crossed pair QUALITY score');

subplot(2,2,4);plot(crossed_pair_scores_twist, crossed_pair_quality_scores_twist,'.')
xlabel( 'Crossed pair score');
ylabel( 'Crossed pair QUALITY score');
