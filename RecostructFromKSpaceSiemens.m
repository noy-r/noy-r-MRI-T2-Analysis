function rec_img = RecostructFromKSpaceSiemens(k_data) 
    fft_dims = [1 3];
        for f = fft_dims
            k_data = ifftshift(ifft(fftshift(k_data, f), [], f), f);
        end
    
        % Combine coils using sum-of-squares (SOS)
        rec_img = squeeze(sqrt(sum(abs(k_data).^2, 2)));
end

