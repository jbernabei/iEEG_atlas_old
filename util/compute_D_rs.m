function D_rs = compute_D_rs(scores,res_elec_inds)

    res_result = scores(res_elec_inds);
    
    non_res_result = scores;
    non_res_result(res_elec_inds) = [];

    mwu_result = mwwtest(res_result',non_res_result');

    D_rs = max(mwu_result.U(2))./(length(res_result).*length(non_res_result));

end