%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cub_info : Load data cub (format .pds)
%
%Description:
%	Load a description file ( .lbl) of a cub
%       give a struct usable by other cub_* functions like cub_read
%
%See also:
%	cub_read cub_view
%
% version 1.05
%Author: Etienne de Foras - Laboratoire de Planetologie de Grenoble
function cubinfo=cub_info(cubename)

  %  [info_file,err_file,msg_file]=stat(cubename);
  %  if err_file ~= 0
%	disp('Error:');
%	disp(cubename);
%	disp(msg_file);
%	cubinfo=-1;
%	return;
 %   end


	%lit juste le header
	[key,value]=read_keys(cubename,0);
	
	record_type=cub_get_field(key,value,'RECORD_TYPE');
	if index(record_type,'FIXED_LENGTH')==0
		error('Unable to handle this RECORD_TYPE !!! Please fixme !');
	end

	record_bytes=cub_get_field(key,value,'RECORD_BYTES');
	label_record=cub_get_field(key,value,'LABEL_RECORDS');
	header_size=str2num(record_bytes)*str2num(label_record); % pour ne pas lire trop loin !

	%lit tout le header, pas plus loin
	[key,value]=read_keys(cubename,header_size);
	cubinfo.key=key;
	cubinfo.value=value;
	
	%lit la taille et le kind
	cubinfo.core.type=cub_get_field(key,value,'CORE_ITEM_TYPE');
	cubinfo.core.elem_size=str2num(cub_get_field(key,value,'CORE_ITEM_BYTES'));
	cubinfo.core.kind=cub_get_field(key,value,'AXIS_NAME');

	if index(cubinfo.core.kind,'(SAMPLE,LINE,BAND)') ~= 0 %BSQ
		v=cub_get_field(key,value,'CORE_ITEMS');
		tmp=sscanf(v,'(%f,%f,%f)',3);
		cubinfo.core.sample=tmp(1);
		cubinfo.core.line=tmp(2);
		cubinfo.core.band=tmp(3);
		v=cub_get_field(key,value,'SUFFIX_ITEMS');
		tmp=sscanf(v,'(%f,%f,%f)',3);
		cubinfo.suffix.sample=tmp(1);
		cubinfo.suffix.line=tmp(2);
		cubinfo.suffix.band=tmp(3);
	end

	if index(cubinfo.core.kind,'(BAND,SAMPLE,LINE)') ~= 0 %BIP
		v=cub_get_field(key,value,'CORE_ITEMS');
		tmp=sscanf(v,'(%f,%f,%f)',3);
		cubinfo.core.band=tmp(1);
		cubinfo.core.sample=tmp(2);
		cubinfo.core.line=tmp(3);
		v=cub_get_field(key,value,'SUFFIX_ITEMS');
		tmp=sscanf(v,'(%f,%f,%f)',3);
		cubinfo.suffix.band=tmp(1);
		cubinfo.suffix.sample=tmp(2);
		cubinfo.suffix.line=tmp(3);
	end
	
	if index(cubinfo.core.kind,'(SAMPLE,BAND,LINE)') ~= 0 %BIL
		v=cub_get_field(key,value,'CORE_ITEMS');
		tmp=sscanf(v,'(%f,%f,%f)',3);
		cubinfo.core.sample=tmp(1);
		cubinfo.core.line=tmp(3);
		cubinfo.core.band=tmp(2);
		v=cub_get_field(key,value,'SUFFIX_ITEMS');
		tmp=sscanf(v,'(%f,%f,%f)',3);
		cubinfo.suffix.sample=tmp(1);
		cubinfo.suffix.line=tmp(3);
		cubinfo.suffix.band=tmp(2);
	end

	%lit les lambda
	slambda=(cub_get_field(key,value,'BAND_BIN_CENTER'));
     
     	if  ~isempty(slambda)
		cubinfo.lambda=[];
    		slambda=slambda(2:end); % pour sauter le '('
		for ihl=1:cubinfo.core.band
			[lam,n]=sscanf(slambda,'%f',1);
        		cubinfo.lambda=[cubinfo.lambda;lam];
			index_virgule=index(slambda,',');
			slambda=slambda(index_virgule+1:end);
		end
	end
	
	%lit le nom de fichier de data
	v=cub_get_field(key,value,'^QUBE');
	if index(v,'"') ~= 0 
		%donnee separee
		simplename=v(index(v,'"')+1:rindex(v,'"')-1);
		last_slash=rindex(cubename,'/');
		if last_slash ~= 0
			simplename=[cubename(1:last_slash),simplename];
		end
		cubinfo.cubename=simplename;
		cubinfo.decall=0;
	else
		%donnee incluse
		cubinfo.decall=(str2num(v)-1)*str2num(record_bytes);
		cubinfo.cubename=cubename;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val=cub_get_field(key,value,field)

    % to replace the cellidx function, not always present
    pos=0;
    %for i=1:columns(key)
    for i=1:numel(key)
	if strcmp(key{i},field)
	    pos=i;
	    break;
	end
    end
%    [pos,err]=cellidx(key,field);
    
    if pos ~= 0
	val=value{pos};
    else
	val=[];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [key,value]=read_keys(filename,read_to_max)
	virgule=',';
	key={};
	value={};
	bRecordType=0;
	bRecordBytes=0;
	bLabelRecords=0;
	
	%lit le .lbl pour extraire les info
	lbl=fopen(filename,'rb');
	
	if(lbl==-1)
		error(['cub_info : Inexistent file: ',filename]);
	end
	
	%pour toutes les lignes lit les champs key=value
	line=fgetl(lbl);
	while ischar(line)

		%enleve les espaces
    	line=line(line ~= ' ');
		pos_egal=index(line,'=');
		if pos_egal ~= 0
			local_key=line(1:pos_egal-1);
			key{end+1}=deblank(local_key); % a verifier l'importance du deblank
			local_value=line(pos_egal+1:end);
			
			if index(local_key,'RECORD_TYPE')~=0
				bRecordType=1;
			end
			if index(local_key,'RECORD_BYTES')~=0
				bRecordBytes=1;
			end
			if index(local_key,'LABEL_RECORDS')~=0
				bLabelRecords=1;
			end
			
			%regarde les lignes multiples
			pos_parenthese_ouvrante=index(line,'(');
			if pos_parenthese_ouvrante ~= 0
				pos_parenthese_fermante=index(line,')');
				while pos_parenthese_fermante==0
						line=fgetl(lbl); %on suppose !=eol
						line=line(line~=' ');
						local_value=strcat(local_value,line);
						pos_parenthese_fermante=index(line,')');
				end
			end

			value{end+1}=deblank(local_value); %a verifier l'importance de deblank
			
			if(read_to_max==0) && (bRecordType*bRecordBytes*bLabelRecords==1);
				break
			end
	
        else
            % pas de egal
		end

		%verifie que l'on est pas a la fin du header
		if(ftell(lbl)>=read_to_max) && read_to_max>0
			break; %fin du header
		end
	
		line=fgetl(lbl);

	end

	fclose(lbl);
end

%%%%%%%% Wrappers of strfind into index and rindex %%%%%
function k = index(str, pattern)
    k = strfind(str, pattern);
    if~isempty(k)
        k=k(1);
    else
        k=0;
    end
end

function k = rindex(str, pattern)
    k = strfind(str, pattern);
    if~isempty(k)
        k=k(numel(k));
    else
        k=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
