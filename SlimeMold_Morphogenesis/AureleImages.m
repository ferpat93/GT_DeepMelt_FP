% Testing Aurele's images
disp_folder='C:\Users\lfp3\Downloads\AI\Image_analysis'; %% Folder where you have the images
srcFiles = dir([disp_folder '/' '*.cr2']);
nColors = 2; %SM - border - background

for i=1:numel(srcFiles)
    rgb_im=imread(fullfile(disp_folder,srcFiles(i).name));
    rgb_im=rgb_im(1000:2600,1000:4500,:);
    
    im=rgb2lab(rgb_im);
    im=im2single(im(:,:,2:3));

    KM_ab = imsegkmeans(im,nColors,'NumAttempts',3);
    KM_rgb = imsegkmeans(rgb_im,nColors,'NumAttempts',3);
    
    figure('Renderer', 'painters', 'Position', [10 10 500 750])
    subplot(3,1,1)
        imshow(rgb_im)
        title('Original Image')
    subplot(3,1,2)
        imagesc(KM_ab)
        colormap(prism(3))
        title('LAB Segmentation')
    subplot(3,1,3)
        imagesc(KM_rgb)
        colormap(prism(3))
        title('RGB Segmentation')
        
    sgtitle(srcFiles(i).name, 'Interpreter', 'none')    
    saveas(gcf,fullfile(disp_folder,strcat(srcFiles(i).name,'-Segmentation.tif')))   
end