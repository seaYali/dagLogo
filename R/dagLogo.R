#' Visualize daglogo using a heatmap
#' 
#' Using a heat map to visualize results of testing differential amino acid uasage
#' 
#' @param testDAUresults An object of \code{testDAUresults} (results of testing 
#' differential amino acid uasage).
#' @param type A character vector of length 1. Type of statistics to display 
#' on y-axis, available choices are "diff" or "zscore".
#' @param ... Other parameters passed to \code{pheatmap} function.
#' @import pheatmap
#' @return The output from \code{pheatmap}
#' @export
#' @author Jianhong Ou
#' @examples
#' data("seq.example")
#' data("proteome.example")
#' bg <- buildBackgroundModel(seq.example, proteome=proteome.example, permutationSize=10)
#' t0 <- testDAU(seq.example, bg)
#' dagHeatmap(testDAUresults = t0, type = "diff")

dagHeatmap <-function(testDAUresults, type = c("diff", "zscore"), ...) 
{
    if (missing(testDAUresults) || class(testDAUresults) != "testDAUresults") 
    {
        stop("testDAUresults should be an object of testDAUresults\n
             Please try ?testDAU to get help.", call. = FALSE)
    }
    type <- match.arg(type)
    dat <- ifelse(type == "diff", testDAUresults@difference, testDAUresults@zscore )
    dat <- dat[order(rowSums(dat), decreasing = TRUE), ]
    pheatmap(dat,
             ...,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             scale = "column")
}

#' Color sets for visualization
#' 
#' Create color encoding for visualization of a peptide sequence logo.
#'
#' @param colorScheme A character vecto of length 1, determining the color scheme
#' based on amino acid classification methods. Available color schemes are "null",
#' "classic", "charge", "chemistry", and "hydrophobicity".
#'
#' @return A named character vector of colors
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' colorsets("classic")

colorsets <-function(colorScheme = c("null", "classic", "charge", "chemistry", 
                                     "hydrophobicity")) 
{
    colorScheme <- match.arg(colorScheme)

    switch(colorScheme,
           null = get("auto", envir = cachedEnv),
           classic = get("classic", envir = cachedEnv)$color,
           charge = get("charge", envir = cachedEnv)$color,
           chemistry = get("chemistry", envir = cachedEnv)$color,
           hydrophobicity = get("hydrophobicity", envir = cachedEnv)$color)
}

#' Get character symbols for grouped amino acids
#'
#' @param nameScheme A character vecto of length 1, determining the character
#' symbols used to represent amino acids grouped by their physical and chemical
#' properties. The available \code{nameScheme} are "classic", "charge", 
#' "chemistry", and "hydrophobicity".
#'
#' @return A named character vector of character symbols
#' @export
#' @author Jianhong Ou, Haibo Liu
#'
#' @examples
#' getNameHash("charge")

getNameHash <-function(nameScheme = c("null", "classic", "charge", "chemistry", 
                                      "hydrophobicity"))
{
        nameScheme <- match.arg(nameScheme)
        switch(
            nameScheme,
            null = NULL,
            classic = get("classic", envir = cachedEnv)$symbol,
            charge = get("charge", envir = cachedEnv)$symbol,
            chemistry = get("chemistry", envir = cachedEnv)$symbol,
            hydrophobicity = get("hydrophobicity", envir = cachedEnv)$symbol
        )
}

#' Create sequence logo
#' 
#' Create sequence logo for visualizing results of testing differential usage 
#' of amino acids.
#'
#' @param testDAUresults An object of \code{testDAUresults} (results of testing 
#' differential amino acid uasage).
#' @param type A character vector of length 1. Type of statistics to display 
#' on y-axis, available choices are "diff" or "zscore".
#' @param pvalueCutoff A numeric vector of length 1. A cutoff of p-value.
#' @param namehash A named character vector, with 3-letter symbols of amino acids
#' as names and single-letter symbols as values.
#' @param font A character vector of length 1. Font type for displaying sequence
#' Logo
#' @param legend A logical vector of length 1, indicating whether to show the 
#' legend.
#' @param labelRelativeToAnchor A logical vector of length 1, indicating whether
#' x-axis label should be adjusted relative to the anchoring position.
#' @param labels A character vector, x-axis labels
#' @param gscmd  A character vector of length 1, indicating the file path to 
#' Ghostscript. In MAC OS, it usually is "/usr/local/bin/gs". This should be 
#' defined as a system variable or provide directly.
#' @param fontface An integer, fontface of text for axis annotation and legends.
#' @param fontsize An integer, fontsize of text for axis annotation and legends.
#'
#' @return 
#' @export
#' @author Jianhong Ou, Haibo Liu
#'
#' @examples
#' data("seq.example")
#' data("proteome.example")
#' bg <- buildBackgroundModel(seq.example, proteome=proteome.example, permutationSize=10)
#' t0 <- testDAU(seq.example, bg)
#' t1 <- testDAU(seq.example, bg, group="classic")
#' t2 <- testDAU(seq.example, bg, group="charge")
#' t3 <- testDAU(seq.example, bg, group="chemistry")
#' t4 <- testDAU(seq.example, bg, group="hydrophobicity")
#' dagLogo(t0)
#' dagLogo(t1, namehash = getNameHash(t1@group))
#' dagLogo(t2, namehash = getNameHash(t2@group))
#' dagLogo(t3, namehash = getNameHash(t3@group))
#' dagLogo(t4, namehash = getNameHash(t4@group))


dagLogo <- function(testDAUresults,
                    type = c("diff", "zscore"),
                    pvalueCutoff = 0.05,
                    namehash = getNameHash(testDAUresults@group),
                    font = "Helvetica-Bold",
                    fontface = "bold",
                    fontsize= 5,
                    legend = FALSE,
                    labelRelativeToAnchor = FALSE,
                    gscmd = Sys.getenv("R_GSCMD"),
                    labels = NULL) 
{

    if (missing(testDAUresults) || class(testDAUresults) != "testDAUresults") 
    {
        stop("testDAUresults should be an object of testDAUresults\n
            Please try ?testDAU to get help.", call. = FALSE)
    }
    type <- match.arg(type)
    if (type == "diff") 
    {
        dat <- testDAUresults@difference
    } else
    {
        dat <- testDAUresults@zscore
    }
    
    npos <- ncol(dat)
    ncha <- nrow(dat)
    if (!is.null(labels)) 
    {
        if (length(labels) < npos) 
        {
            stop("The length of labels is too short!")
        }
    }
    colset <- colorsets(testDAUresults@group)
    colset <- colset[rownames(dat)]
    if (any(is.na(colset)))
        stop("Not every symbol has its color setting.",call. = FALSE)
    if (!is.null(namehash)) 
    {
        rownames(dat) <- namehash[rownames(dat)]
        names(colset) <- namehash[names(colset)]
    }
    rname <- rownames(dat)
    if (max(nchar(rname)) > 1) 
    {
        stop("Please using the namehash to convert the symbols into single letters.",
             call. = FALSE)
    }
    key <- paste("x",
                 ncha,
                 font,
                 paste(colset[rname], collapse = ""),
                 paste(rname, collapse = ""),
                 sep = "_")
    symbolsCache <- ifelse(exists("tmp_motifStack_symbolsCache", envir = cachedEnv), 
                get("tmp_motifStack_symbolsCache", envir = cachedEnv), list())

    if (!is.null(symbolsCache[[key]])) 
    {
        symbols <- symbolsCache[[key]]
    } else 
    {
        symbols <- motifStack:::coloredSymbols(ncha, font, colset[rname], rname)
        symbolsCache[[key]] <- symbols
        assign("tmp_motifStack_symbolsCache", symbolsCache, envir = cachedEnv)
    }
    
    dw <- ifelse(legend, 1 / (npos + 6), 1 / (npos + 2))
    
    ##check font height to fix the number size
    ## plot.new()  ## This doesn't refresh the setting for the new page of 
    ##                the grid graphics system
    grid.newpage()  ## update the plotting frame for the grid graphics system
    line1 <-
        as.numeric(convertUnit(unit(1, "strwidth", "W"), "npc"))
    cex <- dw / line1
    lwd <- cex / 3
    x0 <- dw
    x1 <- 4 / 5 * dw
    x2 <- 6 / 5 * dw
    x3 <- 1 / 5 * dw
    ##set ylim
    ## below x-axis
    datN <- apply(dat, 2, function(.col)
        sum(.col[.col < 0]))
    
    ## above x-axis
    datP <- apply(dat, 2, function(.col)
        sum(.col[.col > 0]))
    
    ylim <-
        c((as.integer(min(datN) / 0.05) - 1) * 0.05, 
          (as.integer(max(datP) / 0.05) + 1) * 0.05)
    
    remap <- function(x) {
        (ylim[1] - x) / (ylim[1] - ylim[2]) / (1 + dw)
    }
    reheight <- function(x) {
        abs(x / (ylim[1] - ylim[2]) / (1 + dw))
    }
    
    ##draw axis
    x.pos <- 0
    grid.text(0,
              x0,
              remap(0) + dw / 2,
              just = c(.5, .5),
              gp = gpar(fontsize=fontsize, fontface = fontface)
    )
    grid.lines(
        x = x0,
        y = c(remap(0) + dw, 1),
        arrow = arrow(length = unit(0.01, "npc")),
        gp = gpar(lwd = lwd/2)
    )
    grid.lines(
        x = x0,
        y = c(remap(0), 0),
        arrow = arrow(length = unit(0.01, "npc")),
        gp = gpar(lwd = lwd/2)
    )
    
    ##draw tick
    tick <- ifelse(type == "diff", 0.1, 10 ^ as.integer(log10(max(
            abs(as.numeric(dat))))))
    
    times <- ifelse(type == "diff", 100, 1)
    for (i in c(as.integer(min(datN) / tick):(-1), 1:as.integer(max(datP) / tick))) 
    {
        grid.lines(
            x = c(x2, x0),
            y = remap(i * tick) + dw / 2,
            gp = gpar(lwd = lwd/2)
        )
        grid.text(
            label = times * tick * i,
            x = x1,
            y = remap(i * tick) + dw / 2,
            just = c(1, .5),
            gp = gpar(fontsize=fontsize, fontface = fontface)
        )
    }
    
    ## add ylab
    grid.text(label= ifelse(type == "diff", "DAU (%)", "Z-score"), 
              x = x3, y = remap(0) + dw / 2,
              just = "centre",  rot = 90, 
              gp = gpar(fontsize=fontsize, fontface = fontface))
    
    x.pos <- x0 + dw / 2
    for (j in 1:npos) {
        heights <- dat[, j]
        id <- order(heights)
        heights <- heights[testDAUresults@pvalue[, j] < pvalueCutoff]
        id <-
            id[id %in% which(testDAUresults@pvalue[, j] < pvalueCutoff)]
        id1 <- order(heights)
        y.pos <- remap(sum(heights[heights < 0]))
        flag <- 0
        if (length(heights) > 0) {
            for (i in 1:length(heights)) {
                h <- reheight(heights[id1[i]])
                if (heights[id1[i]] > 0)
                    flag <- flag + 1
                if (flag == 1) 
                {
                    ## plot x-axis tick labels
                    if (labelRelativeToAnchor) 
                    {
                        grid.text(
                            j - 1 - testDAUresults@upstream,
                            x.pos + dw / 2,
                            y.pos + dw / 2,
                            just = c(.5, .5),
                            gp = gpar(fontsize=fontsize, fontface = fontface)
                        )
                    } else
                    {
                        
                        if (!is.null(labels)) 
                        {
                            grid.text(
                                labels[j],
                                x.pos + dw / 2,
                                y.pos + dw / 2,
                                just = c(.5, .5),
                                gp = gpar(fontsize=fontsize, fontface = fontface)
                            )
                        } else
                        {
                            grid.text(
                                j,
                                x.pos + dw / 2,
                                y.pos + dw / 2,
                                just = c(.5, .5),
                                gp = gpar(fontsize=fontsize, fontface = fontface)
                            )
                        }
                    }
                    y.pos <- y.pos + dw
                    flag <- flag + 1
                }
                ## plot symbols for amino acids 
                if (h > 0) {
                    grid.draw(
                        grImport::pictureGrob(
                            symbols[[id[i]]],
                            x.pos,
                            y.pos,
                            dw,
                            h,
                            just = c(0, 0),
                            distort = TRUE
                        )
                    )
                    y.pos <- y.pos + h
                }
            }
            
        }
        ## plot x-axis tick labels
        if (flag == 0)
        {
            grid.text(j,
                      x.pos + dw / 2,
                      y.pos + dw / 2,
                      just = c(.5, .5),
                      gp = gpar(fontsize=fontsize, fontface = fontface))
        }    
        x.pos <- x.pos + dw
    }
    ## plot legend
    if (legend) 
    {
        if (is.null(namehash)) 
        {
            namehash <- get("namehash", envir = cachedEnv)
        }
        for (i in 1:length(namehash)) 
        {
            grid.text(
                namehash[i],
                x = (npos+2)*dw,
                y = .95 - i * 0.05 * lwd/2,
                just = c(.5, .5),
                gp = gpar(col = colset[namehash[i]], fontsize=fontsize, fontface = fontface))
            grid.text(
                names(namehash)[i],
                x = (npos+2)*dw + 0.03 * lwd/2,
                y = .95 - i * 0.05 * lwd/2,
                just = c(0, .5),
                gp = gpar(col = colset[namehash[i]], fontsize=fontsize, fontface = fontface))
        }
    }
}