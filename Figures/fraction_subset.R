# Fraction calculation from Ludwig

absdiff <- abs(x - y)
ind <- order(absdiff)
fract <- 0.8   # the top fraction that you want to consider
length.out <- round(fract * length(absdiff))
fract.grid <- seq_len(length.out)
plot(x[fract.grid], y[fract.grid])

# Example ----------------------------------------------------------------------
# Assuming your data vectors are both of length 100
#
# absdiff <- abs(x - y)
# ind <- order(absdiff)
# plot(x[ind[1:80]], y[ind[1:80]])  