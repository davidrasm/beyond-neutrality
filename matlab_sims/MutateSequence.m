function output_seq = MutateSequence(input_seq, n_mutations)

% assume HKY model of sequence evolution
transition_transversion_ratio = 3;
%transition_transversion_ratio = 3, so kappa = 2*transition_transversion_ratio = 6.

output_seq = input_seq;

seq_length = length(input_seq);
for i = 1:n_mutations
    rand_site = randi(seq_length);
    this_nuc = output_seq(rand_site);
    if isequal(this_nuc, 'a')
        if rand < transition_transversion_ratio/(transition_transversion_ratio + 1) % then transition
            output_seq(rand_site) = 'g'; 
        else  % then transversion
            if rand < 0.5
                output_seq(rand_site) = 'c'; 
            else
                output_seq(rand_site) = 't'; 
            end
        end
    elseif isequal(this_nuc, 'g')
        if rand < transition_transversion_ratio/(transition_transversion_ratio + 1) % then transition
            output_seq(rand_site) = 'a'; 
        else  % then transversion
            if rand < 0.5
                output_seq(rand_site) = 'c'; 
            else
                output_seq(rand_site) = 't'; 
            end
        end
    elseif isequal(this_nuc, 'c')
        if rand < transition_transversion_ratio/(transition_transversion_ratio + 1) % then transition
            output_seq(rand_site) = 't'; 
        else  % then transversion
            if rand < 0.5
                output_seq(rand_site) = 'a'; 
            else
                output_seq(rand_site) = 'g'; 
            end
        end
    elseif isequal(this_nuc, 't')
        if rand < transition_transversion_ratio/(transition_transversion_ratio + 1) % then transition
            output_seq(rand_site) = 'c'; 
        else  % then transversion
            if rand < 0.5
                output_seq(rand_site) = 'a'; 
            else
                output_seq(rand_site) = 'g'; 
            end
        end
    else
        error('nucleotide has to be a, c, g, or t');
    end
end
