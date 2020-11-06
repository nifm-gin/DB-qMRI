function [] = FormatDico(pre_filename, post_filename, dico_filename)


Pre     = loadjson(pre_filename);
Post    = loadjson(post_filename);

Dico.MRSignals      = abs(Post.MRSignals) ./ abs(Pre.MRSignals); % Ratio post/pre signals 
% Dico.MRSignals      = [abs(Pre.MRSignals) abs(Post.MRSignals)];
Dico.Tacq           = Pre.Sequence.Tacq;
Dico.Parameters.Par = Pre.Parameters.Par; % Parameters used to simulate X signals
Dico.Parameters.Labels = Pre.Parameters.Labels;
clear Pre Post

save(dico_filename,'Dico')