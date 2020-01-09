%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cub_read : Read a slice of a cube (format .pds)
%
%Description:
%	Read a slice of a cub using a cubinfo source
%
% plane=cub_read(cubinfo,linemin,linemax)
% if just cubinfo is given, read all the scanlines
% give back the hyperspectral image as a ND matrix
% 
%See also:
%	cub_info cub_view
%
% version 1.07
%Author: Etienne de Foras - Laboratoire de Planetologie de Grenoble
function cub=cub_read2(cubinfo,linemin,linemax)
	
	% si un seul argument, lit toutes les lignes
	if nargin==1
		linemin=1;
		linemax=cubinfo.core.line;
	end

	elem_size=-1;
	
	if (index(cubinfo.core.type,'PC_REAL')~=0) && (cubinfo.core.elem_size==4 )
		elem_size=4;
		elem_name='real*4';
	end
		
	if (index(cubinfo.core.type,'real*4')~=0) && (cubinfo.core.elem_size==4 )
		elem_size=4;
		elem_name='real*4';
	end

	if (index(cubinfo.core.type,'PC_INTEGER')~=0) && (cubinfo.core.elem_size==2 ) 
		elem_size=2;
		elem_name='integer*2';
	end
	
	if (index(cubinfo.core.type,'LSB_INTEGER')~=0) && (cubinfo.core.elem_size==2 )
		elem_size=2;
		elem_name='integer*2';
	end
	
	if (index(cubinfo.core.type,'LSB_INTEGER')~=0) && (cubinfo.core.elem_size==4 )
		elem_size=4;
		elem_name='integer*4';
	end
	
	if (index(cubinfo.core.type,'LSB_INTEGER')~=0) && (cubinfo.core.elem_size==1 )
		elem_size=1;
		elem_name='integer*1';
	end
	if (index(cubinfo.core.type,'uint8')~=0) && (cubinfo.core.elem_size==1 )
		elem_size=1;
		elem_name='integer*1';
	end
	if elem_size==-1
		error(['cub_read : elem size not supported yet : ',cubinfo.core.type]);
	end

	fp=fopen(cubinfo.cubename,'rb');
	if fp==-1
		error(['cub_read : Unable to open file, !!!!',cubinfo.cubename]);
    end
	
	if strcmp(cubinfo.core.kind,'(SAMPLE,LINE,BAND)') %BSQ
		fseek(fp,cubinfo.decall,'bof');
		cub2=fread(fp,cubinfo.core.sample*cubinfo.core.line*(cubinfo.core.band+cubinfo.suffix.band),elem_name,0);
		cub2=reshape(cub2,[cubinfo.core.sample,cubinfo.core.line,cubinfo.core.band+cubinfo.suffix.band]);	
		if ndims(cub2)==3
			cub2=permute(cub2,[2,1,3]);
		else
			cub2=permute(cub2,[2,1]);
		end
	end

	if strcmp(cubinfo.core.kind,'(BAND,SAMPLE,LINE)')==1 %BIP
		fseek(fp,cubinfo.decall+elem_size*cubinfo.core.sample*(linemin-1)*(cubinfo.core.band+cubinfo.suffix.band),'bof');
		cub2=fread(fp,[cubinfo.core.sample*(linemax-linemin+1)*(cubinfo.core.band+cubinfo.suffix.band)],elem_name,0);
		cub2=reshape(cub,[cubinfo.core.band+cubinfo.suffix.band,cubinfo.core.sample,(linemax-linemin+1)]);
		cub2=permute(cub2,[3,2,1]);
	else

	if strcmp(cubinfo.core.kind,'(SAMPLE,BAND,LINE)') %BIL
		fseek(fp,cubinfo.decall,'bof');
		cub2=fread(fp,cubinfo.core.sample*cubinfo.core.line*(cubinfo.core.band+cubinfo.suffix.band),elem_name,0);
		cub2=reshape(cub,[cubinfo.core.sample,cubinfo.core.band+cubinfo.suffix.band,cubinfo.core.line]);
		cub2=permute(cub2,[3,1,2]);
	end
		%else ? -> error ?
	end

	fclose(fp);
	
	%enleve les suffixes et les suffixes
	cub2=cub2(linemin:linemax,1:cubinfo.core.sample,1:cubinfo.core.band);	
    cub=struct;
	cub.spectels=cub2;
	cub.spectels_flat=[];%	cub_flatten(cub.spectels);
	cub_size=size(cub.spectels);
	cub.line=cub_size(1);
	cub.sample=cub_size(2);
	if ndims(cub)==3
		cub.band=cub_size(3);
	else
		cub.band=1;
	end
	cub.masked=false;
%	cub=cub_update(cub);
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

% function k = rindex(str, pattern)
%     k = strfind(str, pattern);
%     if~isempty(k)
%         k=k(numel(k));
%     else
%         k=0;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%