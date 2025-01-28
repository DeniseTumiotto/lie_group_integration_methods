function rslt = qpvkp(p,v)

rslt = qIm(qp(qp(p, [0;v]), qconj(p)));