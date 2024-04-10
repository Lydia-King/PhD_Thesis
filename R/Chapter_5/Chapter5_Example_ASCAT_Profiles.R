# Chapter 5: ASCAT Profiles Example

## Single Plot functions
base.gw.plot = function(nAfullPlot,
                        nBfullPlot,
                        colourTotal,
                        colourMinor,
                        y_limit,
                        twoColours = FALSE,
                        main) {
    par(
        mar = c(2.5, 5, 5, 0.5),
        cex = 0.4,
        cex.main = 3,
        cex.axis = 2.5
    )
    ticks = seq(0, y_limit, 1)
    plot(
        c(1, length(
            rep(nAfullPlot$values, nAfullPlot$lengths)
        )),
        c(-0.06, y_limit),
        type = "n",
        xaxt = "n",
        yaxt = "n",
        main = main,
        xlab = "",
        ylab = ""
    )
    axis(side = 2, at = ticks)
    abline(h = ticks,
           col = "lightgrey",
           lty = 1)
    axis(side = 1, at = seq(0, length(
        rep(nAfullPlot$values, nAfullPlot$lengths)
    ), 50000000))
    
    
    A_rle <- nAfullPlot
    start = 0
    # plot total copy number
    for (i in 1:length(A_rle$lengths)) {
        val <- A_rle$values[i]
        size <- A_rle$lengths[i]
        rect(
            start,
            (val - 0.072),
            (start + size - 1),
            (val + 0.072),
            col = ifelse((twoColours &
                              val >= y_limit),
                         adjustcolor(
                             colourTotal,
                             red.f = 0.75,
                             green.f = 0.75,
                             blue.f = 0.75
                         ),
                         colourTotal
            ),
            border = ifelse((twoColours &
                                 val >= y_limit),
                            adjustcolor(
                                colourTotal,
                                red.f = 0.75,
                                green.f = 0.75,
                                blue.f = 0.75
                            ),
                            colourTotal
            )
        )
        start = start + size
    }
    
    B_rle <- nBfullPlot
    start = 0
    # plot minor allele copy number
    for (i in 1:length(B_rle$values)) {
        val <- B_rle$values[i]
        size <- B_rle$lengths[i]
        rect(
            start,
            (val - 0.072),
            (start + size - 1),
            (val + 0.072),
            col = ifelse((twoColours &
                              val >= y_limit),
                         adjustcolor(
                             colourMinor,
                             red.f = 0.75,
                             green.f = 0.75,
                             blue.f = 0.75
                         ),
                         colourMinor
            ),
            border = ifelse((twoColours &
                                 val >= y_limit),
                            adjustcolor(
                                colourMinor,
                                red.f = 0.75,
                                green.f = 0.75,
                                blue.f = 0.75
                            ),
                            colourMinor
            )
        )
        start = start + size
    }
    
    chrk_tot_len = sum(B_rle$lengths)
    tpos = (chrk_tot_len) / 2
    
    abline(v = 0,
           lty = 1,
           col = "lightgrey")
    
}


ascat.plotAscatProfile <- function(n1all, n2all, y_limit = 2.06, main) {
    nA2 = inverse.rle(n1all)
    nB2 = inverse.rle(n2all)
    nA <- nA2
    nB <- nB2
    
    nBPlot <- ifelse(nB <= y_limit, nB + 0.072, y_limit + 0.072)
    nAPlot <- ifelse(nA <= y_limit, nA - 0.072, y_limit + 0.072)
    
    colourTotal = "#E03546" # red
    colourMinor = "#3557E0" # blue
    
    maintitle = ""
    base.gw.plot(
        rle(nAPlot),
        rle(nBPlot),
        colourTotal,
        colourMinor,
        main,
        y_limit = y_limit,
        twoColours = TRUE
    )
}

## Example Profiles
dataset1_B <-
    list(values = c(1), lengths = c(250000000)) 

dataset1_A <-
    list(values = c(1, 2, 1, 1),
         lengths = c(50000000, 17208303, 20000000, (250000000 - (
             50000000 + 17208303 + 20000000
         ))))

dataset1_A_2 <-
    list(
        values = c(1, 1, 2, 0, 1),
        lengths = c(50000000, 100000000, 22105618, 11299836, (
            250000000 - (50000000 + 100000000 + 22105618 + 11299836)
        ))
    )

dataset1_C <-
    list(
        values = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
        lengths = c(
            20000000,
            12466503,
            10000000,
            12466503,
            15000000,
            12466503,
            13000000,
            12466503,
            11299836,
            12466503,
            (250000000 - (sum(
                c(
                    20000000,
                    12466503,
                    10000000,
                    12466503,
                    15000000,
                    12466503,
                    13000000,
                    12466503,
                    11299836,
                    12466503
                )
            )))
        )
    )

dataset1_C1 <-
    list(values = c(1, 2, 0, 1),
         lengths = c(2000000, 22105618 * 3, 11299836 * 3, (
             250000000 - (2000000 + 22105618 * 3 + 11299836 * 3)
         )))

png(filename = "../../figures/Chapter_5/Example_ASCAT_Scen.png",
    width = 900,
    height = 680)
par(mfrow = c(2, 2), oma = c(1, 1, 6, 1))
ascat.plotAscatProfile(dataset1_B, dataset1_B, main = "(A)")
ascat.plotAscatProfile(dataset1_A, dataset1_B, main = "(B)")
ascat.plotAscatProfile(dataset1_A_2, dataset1_B, main = "(C)")
ascat.plotAscatProfile(dataset1_C1, dataset1_C, main = "(D)")
mtext(
    expression(bold(
        "Simulated Allele-specific Copy Number Profiles"
    )),
    line = 0,
    side = 3,
    outer = TRUE,
    cex = 2
)
dev.off()

## Plot for Annotated ASCAT 
ascat.plotAscatProfile1 <- function(n1all, n2all, y_limit = 2.06, main) {
    nA2 = inverse.rle(n1all)
    nB2 = inverse.rle(n2all)
    nA <- nA2
    nB <- nB2
    
    nBPlot <- ifelse(nB <= y_limit, nB + 0.072, y_limit + 0.072)
    nAPlot <- ifelse(nA <= y_limit, nA - 0.072, y_limit + 0.072)
    
    colourTotal = "#E03546" # red
    colourMinor = "white" # blue
    
    maintitle = ""
    base.gw.plot(
        rle(nAPlot),
        rle(nBPlot),
        colourTotal,
        colourMinor,
        main,
        y_limit = y_limit,
        twoColours = TRUE
    )
}

ascat.plotAscatProfile2 <-
    function(n1all, n2all, y_limit = 2.06, main) {
        nA2 = inverse.rle(n1all)
        nB2 = inverse.rle(n2all)
        nA <- nA2
        nB <- nB2
        
        nBPlot <- ifelse(nB <= y_limit, nB + 0.072, y_limit + 0.072)
        nAPlot <- ifelse(nA <= y_limit, nA - 0.072, y_limit + 0.072)
        
        colourTotal = "blue" # red
        colourMinor = "white" # blue
        
        maintitle = ""
        base.gw.plot(
            rle(nAPlot),
            rle(nBPlot),
            colourTotal,
            colourMinor,
            main,
            y_limit = y_limit,
            twoColours = TRUE
        )
    }

dataset1_B <- list(values = c(1), lengths = c(250000000))

dataset1_C1 <- list(
    values = c(1, 2, 0, 1, 2, 1, 0, 1),
    lengths = c(
        20000000,
        22105618 * 3,
        11299836 * 3,
        22105618,
        22105618 * 2,
        35000000,
        2000000,
        (
            250000000 - (
                2000000 + 22105618 * 3 + 11299836 * 3 + 22105618 + 
                    22105618 * 2 + 35000000 + 2000000
            )
        )
    )
)

png(filename = "../../figures/Chapter_5/ASCAT_Profile.png",
    width = 1300,
    height = 480)
par(mfrow = c(1, 2), oma = c(1, 1, 6, 1))
ascat.plotAscatProfile1(dataset1_C1, dataset1_B, main = "Major Allele")
ascat.plotAscatProfile2(dataset1_B, dataset1_B, main = "Minor Allele")
mtext(
    expression(bold(
        "Allele-specific Copy Number Profile Schematic"
    )),
    line = 0,
    side = 3,
    outer = TRUE,
    cex = 2
)
dev.off()
