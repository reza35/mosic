function rgauss, xarr, p
  model = p[0] + p[1]*exp(-((xarr - p[2])/p[3])^2)
  return, model
end
