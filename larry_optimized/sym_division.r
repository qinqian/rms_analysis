df = read.table("barcode_within_tumor_Regular_1.txt")
stat = apply(df, 1, function(x) {
    length(unique(x[!is.na(x)]))
    })

sum(stat == 1)
sum(stat > 1)
