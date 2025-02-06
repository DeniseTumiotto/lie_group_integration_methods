function rslt = T_inv(x)

    norm_x = norm(x);
    f6 = 0;

    if norm_x < 1e-2
        f6 = 1/12 + norm_x^2/720 + norm_x^4/30240;
    else
        f6 = (2 - norm_x * cot(norm_x/2))/2*norm_x^2;
    end

    rslt = eye(3) + skw(x)/2 + f6 * skw(x)^2;

end