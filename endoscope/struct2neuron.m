function neuron_class = struct2neuron(neuron)
%% convert a struct variable to Sources2D class 
% inputs: 
%   neuron: struct variable containing all data values 
% outputs: 
%   neuron_class: Sources2D class variable 
%% author: Pengcheng Zhou, Carnegie Mellon University, 2016 

neuron_class = Sources2D(); 
fields = fieldnames(neuron); 
for m=1:length(fields)
    eval(sprintf('neuron_class.%s=neuron.%s;', fields{m}, fields{m})); 
end