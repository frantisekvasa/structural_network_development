# Formats r-squared and p-values within titles (mainly due to bugs in interaction of rstudio and the b-quote function)

# written by Frantisek Vasa (fv247@cam.ac.uk)

# Pearson's r
rp.main = function(rsq,p,nsig) {
  main = bquote(r^2 ~ '=' ~ .(toString(signif(rsq,nsig))) * ', p =' ~ .(toString(signif(p,nsig)))) 
  return(main)
}

# Spearman's rho
rp.main.sp = function(sp.rho,p,nsig) {
  if (p < 1e-10) {
    main = bquote(rho ~ '=' ~ .(toString(formatC(sp.rho,digits=nsig,format="fg",flag="#"))) * ', p < 1e-10') 
  } else {
    main = bquote(rho ~ '=' ~ .(toString(formatC(sp.rho,digits=nsig,format="fg",flag="#"))) * ', p =' ~ .(toString(signif(p,nsig))))
  }
  return(main)
}