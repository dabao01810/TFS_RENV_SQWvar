function TFS_RENV_savefcn

global def
global work
global set

% Creates or loads .dat file
if exist([def.result_path,work.filename,'.dat'])==2,
    fid = fopen([def.result_path,work.filename,'.dat'],'a');
else
    fid = fopen([def.result_path,work.filename,'.dat'],'w');
	
	fprintf(fid,'Processing: %s\nNoise: %s\nEar: %s\nMode: %s\n',work.userpar3, work.userpar8, work.userpar4, work.userpar5);
	fprintf(fid,'Source Level: %f\nSNR: %f\n',work.userpar7, work.userpar9);
	fprintf(fid,'Loss: %s\nNAL: %i\nSim: %s\n',work.userpar10, work.userpar11, work.userpar12);
	fprintf(fid,'\n\n');
end

percent_correct = 100*sum(work.correct{1})/length(work.correct{1});
confusion = zeros(16,16);
for k=1:length(work.position{1});
	confusion(work.answer{1}(k),work.position{1}(k)) = confusion(work.answer{1}(k),work.position{1}(k)) + 1;
end
confusion = confusion./(ones(size(confusion,1),1)*sum(confusion,1));
confusion = fliplr(flipud(confusion));

% curr_mode = 'Normal';
% switch lower(deblank(work.userpar3))
% 	case {'tfs'},
%         curr_mode = 'TFS';
% 	case {'renv'},
%         curr_mode = 'RENV';
% 	case {'tfs-renv'},
% 		if (rem(work.exppar1,2)==1)
%             curr_mode = 'TFS';
% 		else
%             curr_mode = 'RENV';
% 		end
% 	case {'renv-tfs'},
% 		if (rem(work.exppar1,2)==0)
%             curr_mode = 'TFS';
% 		else
%             curr_mode = 'RENV';
% 		end
% end

fprintf(fid,'Run: %i\n',work.exppar1);
fprintf(fid,'Noise Type: %i\n',work.exppar3);
fprintf(fid, 'SNR: %f\n', work.exppar2);
fprintf(fid,'Processing: %s\n',work.userpar3);
% if ~strcmpi(curr_mode,'RENV')
%     fprintf(fid,'NBands TFS: %i\nNBands RENV: N/A\n', work.exppar2);
% else
%     fprintf(fid,'NBands TFS: %i\nNBands RENV: %i\n', work.exppar2, work.exppar3);
% end
fprintf(fid,'Percent Correct: %4.3f\n',percent_correct);
for k=1:16,
	fprintf(fid, '\t%s',def.consonants{k});
end
fprintf(fid,'\n');
for k1=1:16,
	fprintf(fid,'%s',def.consonants{k1});
	for k2 = 1:16,
		fprintf(fid, '\t%4.3f',confusion(k1,k2));
	end
	fprintf(fid,'\n');
end
fprintf(fid,'\n\n');

fclose(fid);


% Creates or loads .dat file with additional information
if exist([def.result_path,work.filename,'_extra.dat'])==2,
    fid = fopen([def.result_path,work.filename,'_extra.dat'],'a');
else
    fid = fopen([def.result_path,work.filename,'_extra.dat'],'w');
    	
	fprintf(fid,'Processing: %s\nNoise: %s\nEar: %s\nMode: %s\n',work.userpar3, work.userpar8, work.userpar4, work.userpar5);
	fprintf(fid,'Source Level: %f\nSNR: %f\n',work.userpar7, work.userpar9);
	fprintf(fid,'Loss: %s\nNAL: %i\nSim: %s\n',work.userpar10, work.userpar11, work.userpar12);
	fprintf(fid,'\n\n');
end

fprintf(fid,'Run: %i\n',work.exppar1);
fprintf(fid,'Noise Type: %i\n',work.exppar3);
fprintf(fid, 'SNR: %f\n', work.exppar2);
fprintf(fid,'Processing: %s\n',work.userpar3);
fprintf(fid,'\n');

% Add in additional information for VC noise
fprintf(fid, 'Test Stimuli\t Response\t VC Speakers\n');
for k = 1:def.expvarnum{1},
    response = def.consonants{17 - work.answer{1}(k)};
    vc_speakers = 'N/A'
    fprintf(fid, '%s\t\t %s\t\t\t %s\n', work.test_stimuli{k}, response, vc_speakers);
end
fprintf(fid,'\n\n');

fclose(fid);