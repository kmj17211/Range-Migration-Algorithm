function SAR_db_img = SAR_db_Vis(SAR_db_img, X, Y, img_title)
    SAR_db_img = abs(SAR_db_img);
    SAR_db_img = 20 * log10((SAR_db_img - min(SAR_db_img(:))) / (max(SAR_db_img(:)) - min(SAR_db_img(:))));
    SAR_db_img(SAR_db_img < -60) = -60;
    figure; imagesc(Y, X, SAR_db_img); xlabel('Range Distance [m]'); ylabel('Azimuth Distance [m]'); axis xy; title(img_title); colormap('jet'); axis([1.5 4.5 -1 1]); colorbar
end