close all;
clear;
numOfGrids=512;

outDirs={'outA','outB','outC','outD','outE'};

for ii=1:length(outDirs)
    pipimage=[];
    outDir=outDirs{ii};
    
    
    %     load([outDir '/PIP3_PIP2_PTEN.mat']);
    currentFileName=sprintf([outDir '/PIP3.dat']);
    fid=fopen(currentFileName,'r');
    PIP3=fread(fid,inf,'double', 'ieee-le');
    min(min(PIP3))
    MaxCont=max(PIP3);
    PIP3=reshape(PIP3,numOfGrids,length(PIP3)/numOfGrids);
    fclose(fid);
    
    currentFileName=sprintf([outDir '/PTEN.dat']);
    fid=fopen(currentFileName,'r');
    PTEN=fread(fid,inf,'double', 'ieee-le');
    PTEN=reshape(PTEN,numOfGrids,length(PTEN)/numOfGrids);
    fclose(fid);
    
    %%
    maxPTEN=0.9*max(max(PTEN(:,600:end)));
    minPTEN=min(min(PTEN));
    a=1/(maxPTEN-minPTEN);
    b=-a*minPTEN;
    pipimage(:,:,1)=a*PTEN'+b;
    
    maxPIP3=0.9*max(max(PIP3));
    minPIP3=min(min(PIP3));
    a=1/(maxPIP3-minPIP3);
    b=-a*minPIP3;
    pipimage(:,:,2)=a*PIP3'+b;
    
    pipimage(:,:,3)=zeros(size(PIP3'));
    
    %%
    
    
    if ii == 1 || ii == 2 || ii == 3
        %ABC
        t0=451;
    else
        %DE
        t0=1;
    end
    
    %%
    %imwrite(Im16,[outDir '/' outDir '.tif'],'tif','Compression','lzw');
    f1=figure(1);
    
    dt=5;
    disp=t0:1:(t0-1)+512;
    
    imshow(pipimage(disp,1:4:end,:));
    imwrite(pipimage(disp,1:4:end,:),[outDir '/' outDir '.png'],'png');
    %%
end
