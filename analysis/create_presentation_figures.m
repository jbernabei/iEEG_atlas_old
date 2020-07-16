load('iEEG_atlas/output/figure_2A_data.mat')
dlmwrite('data_file.node',good_mean_conn,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','data_file.node','final_render.mat','render_figure.jpg')