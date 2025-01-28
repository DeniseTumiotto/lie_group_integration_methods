function rslt = inv_s3sdr3(x)
    if max(size(x)) == 7
        rslt(1:4) = qconj(x(1:4))/norm(x(1:4));
        rslt(5:7) = 1./(x(5:7));
    else % only s3
        rslt(1:4) = qconj(x(1:4))/norm(x(1:4));
    end

end