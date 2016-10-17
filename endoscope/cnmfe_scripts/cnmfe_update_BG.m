clear Ysignal;
tic; 
Ybg = Y-neuron.A*neuron.C;
rr = ceil(neuron.options.gSiz * bg_neuron_ratio); 
active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
[Ybg, Ybg_weights] = neuron.localBG(Ybg, spatial_ds_factor, rr, active_px, sn, 5); % estiamte local background.
% subtract the background from the raw data.
Ysignal = Y - Ybg;