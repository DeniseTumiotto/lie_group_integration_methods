function rslt = log_s3(q)

    if q(1) > 1 - 5e-3
        q0 = 2 - 2*(q(1)-1)/3 + 4*(q(1)-1)^2/15 - 4*(q(1)-1)^3/35 + 16*(q(1)-1)^4/315 - 16*(q(1)-1)^5/693;
    else
        q0 = 2*acos(q(1))/sqrt(1-q(1)^2);
    end

    rslt = reshape(q0*q(2:4),3,1);
end