function [powerenv] = powerenvelopeConcise(brq,br_ratios,phoneme_sets,phrase_info_num,N)

F_model =  @(x,i,br)(x(1) + x(3).*(x(2)+i)).*exp(x(4).*(x(2) + i)) + br*x(5);

powerenv = zeros(1,N);
for kk =1:length(phrase_info_num(1,:))
        if br_ratios(kk) ~= -1
            br = brq(br_ratios(kk));
        else
            br = 1;
        end
        phoneme_set = phoneme_sets{1,kk};
        for ii=1:length(phoneme_set) 
            xhat =  phoneme_set{1,ii}.xhat;
            n_i = phoneme_set{1,ii}.n_i;
            i =  phoneme_set{1,ii}.i;
            powerenv(n_i)  = F_model(xhat,i,br);

        end
    
end
