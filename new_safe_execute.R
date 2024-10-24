safe_execute <- function(expr, object_name, data_vector) {
  tryCatch(
    {
      # Use bquote to inject data_vector into the expression
      result <- eval(bquote(.(expr)), envir = list(dat = data_vector))
      return(result)
    },
    error = function(e) {
      cat(paste("Error with fitting", object_name, ":", e$message, ";\n",
                "Other methods were still fit.\n"))
      return(NA)
    }
  )
}

testfun <- function(datvec) {
  # Inject the data vector into the Amoroso estimation function
  amo <- safe_execute(quote(estimate_amoroso(dat, plot=0, criterion="maxL")), "amo", datvec)
  return(amo)
}

dat1 <- rgg4(1000, a=4, l=1, c=7, mu=0)
hist(dat1)
testfun(dat1)
