Class 6 R Functions
================
Tuan
1/25/2019

File reading (again!)
---------------------

Here we try to use **read.table()** and friends to unput some example data into R

Click the insert button to insert a code chunk.

``` r
read.table("https://bioboot.github.io/bggn213_S18/class-material/test1.txt", header = T, sep = ",")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
file1 <- "https://bioboot.github.io/bggn213_S18/class-material/test1.txt"
read.csv(file1)
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
file2 <- "https://bioboot.github.io/bggn213_S18/class-material/test2.txt"
read.table(file2, header = T, sep = "$")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
file3 <- "https://bioboot.github.io/bggn213_S18/class-material/test3.txt"
read.table(file3)
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

``` r
add <- function(x, y=1) {
# Sum the input x and y 
  x+y
}

# add(1,2,3). Return error due to unsued argument (the third one)

# add(x=1, y="b"). Return error due to different types of data
```

``` r
rescale <- function(x) {
   rng <-range(x)
   (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
# How would you get your function to work here...
rescale( c(1,2,NA,3,10) )
```

    ## [1] NA NA NA NA NA

``` r
rescale2 <- function(x) {
   rng <-range(x, na.rm = T)
   (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale2(c(1,2, NA, 3, 10))
```

    ## [1] 0.0000000 0.1111111        NA 0.2222222 1.0000000

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
   if(na.rm) {
     rng <-range(x, na.rm=na.rm)
   } else {
     rng <-range(x)
   }
   print("Hello")
   answer <- (x - rng[1]) / (rng[2] - rng[1])
   print("is it me you are looking for?")
   if(plot) {
      plot(answer, typ="b", lwd=4)
     print("please done ever sing again!")
}
   print("I can see it in ...")
}
```

``` r
library(bio3d)
```

``` r
# Can you improve this analysis code?
library(bio3d)
s1 <- read.pdb("4AKE")  # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s2 <- read.pdb("1AKE")  # kinase no drug
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y")  # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")

s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b

plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-13-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-13-3.png)

``` r
hc <- hclust( dist( rbind(s1.b, s2.b, s3.b) ) )
plot(hc)
```

![](class06_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
Bfac.trend <- function(x) {
  library(bio3d) 
  for (i in x) {
    prot <- read.pdb(i)
    prot.chainA <- trim.pdb(prot, chain = "A", elety = "CA")
    prot.b <- prot.chainA$atom$b
    plotb3(prot.b, sse = prot.chainA, type = "l", ylab = "Bfactor")
  }
}
```

``` r
Bfac.trend(c("4AKE", "1AKE", "1E4Y"))
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 0n/dzcbp5b17_l08r6p_l358hsh0000gn/T//RtmpTNa0VU/4AKE.pdb exists. Skipping
    ## download

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 0n/dzcbp5b17_l08r6p_l358hsh0000gn/T//RtmpTNa0VU/1AKE.pdb exists. Skipping
    ## download

![](class06_files/figure-markdown_github/unnamed-chunk-16-1.png)

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## 0n/dzcbp5b17_l08r6p_l358hsh0000gn/T//RtmpTNa0VU/1E4Y.pdb exists. Skipping
    ## download

![](class06_files/figure-markdown_github/unnamed-chunk-16-2.png)![](class06_files/figure-markdown_github/unnamed-chunk-16-3.png)
