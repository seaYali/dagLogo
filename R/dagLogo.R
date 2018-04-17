#' Visualize daglogo using a heatmap.
#' 
#' Using a heat map to visualize results of testing differential amino acid uasage.
#' 
#' @param testDAUresults An object of \code{testDAUresults} (results of testing 
#' differential amino acid uasage).
#' @param type A character vector of length 1, the type of metrics to display 
#' on y-axis. The available options are "diff" and "statistics", which are 
#' differences in amino acide usage at each position between the inputSet and 
#' the backgroundSet, and the Z-scores or odds ratios when Z-test or Fisher's 
#' exact test is performed to test the differential uasge of amino acid at each 
#' position between the two sets.
#' @param ... Other parameters passed to \code{pheatmap} function.
#' @import pheatmap
#' @return The output from \code{pheatmap}
#' @export
#' @author Jianhong Ou
#' @examples
#' data("seq.example")
#' data("proteome.example")
#' bg <- buildBackgroundModel(seq.example, proteome=proteome.example, 
#'                            numSubsamples=10)
#' t0 <- testDAU(seq.example, bg)
#' dagHeatmap(testDAUresults = t0, type = "diff")

dagHeatmap <-function(testDAUresults, type = c("diff", "statistics"), ...) 
{
    if (missing(testDAUresults) || class(testDAUresults) != "testDAUresults") 
    {
        stop("testDAUresults should be an object of testDAUresults\n
             Please try ?testDAU to get help.", call. = FALSE)
    }
    type <- match.arg(type)
    data <- getData(type, testDAUresults)
    dat <- data$dat[order(rowSums(data$dat), decreasing = TRUE), ]
    pheatmap(dat,
             ...,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             scale = "column",
             main = paste0("Heatmap showing ", data$label))
}

#' Get the data for visualization.
#'
#' A helper function to Get the data and it lable for visualization.
#' 
#' @param type A character vector of length 1, the type of metrics to display 
#' on y-axis. The available options are "diff" and "statistics", which are 
#' differences in amino acide usage at each position between the inputSet and 
#' the backgroundSet, and the Z-scores or odds ratios when Z-test or Fisher's 
#' exact test is performed to test the differential uasge of amino acid at each 
#' position between the two sets.
#' @param testDAUresults An object of \code{testDAUresults} (results of testing 
#' differential amino acid uasage).
#'
#' @return A list containing the following compoenets:
#' label A character vector of length 1. The type of data for visualization.
#' dat   A matrix of numeric data for visualization.
#' @export
#' @keywords internal
#' 
#' @author Haibo Liu
#'
#' @examples
#' 
getData <- function(type, testDAUresults)
{
    testType <- testDAUresults@testType
    if (type == "diff")
    {
        dat <- testDAUresults@difference
        label <- "DAU %"
    } else
    {
        dat <- testDAUresults@statistics
        if (testType == "fisher")
        {
            label <- "OddsRatio"
        } else
        {
            label <- "Z-score"
        }
    }
    list(label = label, dat = dat)
}



#' Color sets for visualization.
#' 
#' Create color encoding for visualization of a peptide sequence logo.
#'
#' @param colorScheme A character vecto of length 1, determining the color scheme
#' based on amino acid classification methods. The available \code{colorScheme} are "no",
#' "classic", "charge", "chemistry", "hydrophobicity", "BLOSM50_L1", "BLOSM50_L2", 
#' "BLOSM50_L3", "BLOSM50_L4", "BLOSM50_L5", "BLOSM50_L6", "BLOSM50_L7", 
#' "BLOSM50_L8", "chemistry_property_Mahler", "contact_potential_Maiorov", 
#' "protein_blocks_Rogov", "sequence_alignment_Dayhoff", "structure_alignments_Mirny",
#' and "custom". If "custom" is used, users must define a grouping scheme using a
#' a list containing sublist named as "color", "symbol" and "group" using the 
#' function addGroupingScheme. 
#'
#' @seealso {\code{\link{addGroupingScheme}}}
#' @return A named character vector of colors
#' @export
#' @keywords internal
#' 
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' colorsets("classic")

colorsets <-function(colorScheme = c("no", "classic", "charge", "chemistry", 
                                     "hydrophobicity", "BLOSM50_L1", 
                                     "BLOSM50_L2", "BLOSM50_L3", "BLOSM50_L4", 
                                     "BLOSM50_L5", "BLOSM50_L6", "BLOSM50_L7", 
                                     "BLOSM50_L8", "chemistry_property_Mahler", 
                                     "contact_potential_Maiorov", 
                                     "protein_blocks_Rogov", 
                                     "sequence_alignment_Dayhoff", 
                                     "structure_alignments_Mirny", "custom")) 
{
    colorScheme <- match.arg(colorScheme)

    switch(colorScheme,
           no = get("no", envir = cachedEnv)$color,
           classic = get("classic", envir = cachedEnv)$color,
           charge = get("charge", envir = cachedEnv)$color,
           chemistry = get("chemistry", envir = cachedEnv)$color,
           hydrophobicity = get("hydrophobicity", envir = cachedEnv)$color,
           BLOSM50_L1 = get("BLOSM50_L1", envir = cachedEnv)$color,
           BLOSM50_L2 = get("BLOSM50_L2", envir = cachedEnv)$color,
           BLOSM50_L3 = get("BLOSM50_L3", envir = cachedEnv)$color,
           BLOSM50_L4 = get("BLOSM50_L4", envir = cachedEnv)$color,
           BLOSM50_L5 = get("BLOSM50_L5", envir = cachedEnv)$color,
           BLOSM50_L6 = get("BLOSM50_L6", envir = cachedEnv)$color,
           BLOSM50_L7 = get("BLOSM50_L7", envir = cachedEnv)$color,
           BLOSM50_L8 = get("BLOSM50_L8", envir = cachedEnv)$color,
           chemistry_property_Mahler = get("chemistry_property_Mahler", 
                                           envir = cachedEnv)$color, 
           contact_potential_Maiorov = get("contact_potential_Maiorov",
                                           envir = cachedEnv)$color, 
           protein_blocks_Rogov = get("protein_blocks_Rogov", 
                                      envir = cachedEnv)$color,
           sequence_alignment_Dayhoff = get("sequence_alignment_Dayhoff",
                                            envir = cachedEnv)$color, 
           structure_alignments_Mirny = get("structure_alignments_Mirny", 
                                            envir = cachedEnv)$color,
           custome = get("custom", envir = .globalEnv)$color)
}

#' Get character symbols for grouped amino acids
#'
#' @param groupingScheme A character vecto of length 1, determining the character
#' symbols used to represent amino acids grouped by their physical and chemical
#' properties. The available \code{groupingScheme} are "no",
#' "classic", "charge", "chemistry", "hydrophobicity", "BLOSM50_L1", "BLOSM50_L2", 
#' "BLOSM50_L3", "BLOSM50_L4", "BLOSM50_L5", "BLOSM50_L6", "BLOSM50_L7", 
#' "BLOSM50_L8", "chemistry_property_Mahler", "contact_potential_Maiorov", 
#' "protein_blocks_Rogov", "sequence_alignment_Dayhoff", "structure_alignments_Mirny",
#' and "custom". If "custom" is used, users must define a grouping scheme using a
#' a list containing sublist named as "color", "symbol" and "group" using the 
#' function addGroupingScheme. 
#'
#' @seealso {\code{\link{addGroupingScheme}}}
#'
#' @return A named character vector of character symbols
#' @export
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu
#'
#' @examples
#'

getGroupingSymbol <-function(groupingScheme = c("no", "classic", "charge", 
                                                "chemistry", "hydrophobicity",
                                                "BLOSM50_L1", "BLOSM50_LL2", 
                                                "BLOSM50_L3", "BLOSM50_L4", 
                                                "BLOSM50_L5", "BLOSM50_L6", 
                                                "BLOSM50_L7", "BLOSM50_L8", 
                                                "chemistry_property_Mahler", 
                                                "contact_potential_Maiorov", 
                                                "protein_blocks_Rogov", 
                                                "sequence_alignment_Dayhoff", 
                                                "structure_alignments_Mirny", 
                                                "custom"))
{
    groupingScheme <- match.arg(groupingScheme)
    switch(
        groupingScheme,
        no = NULL,
        classic = get("classic", envir = cachedEnv)$symbol,
        charge = get("charge", envir = cachedEnv)$symbol,
        chemistry = get("chemistry", envir = cachedEnv)$symbol,
        hydrophobicity = get("hydrophobicity", envir = cachedEnv)$symbol,
        BLOSM50_L1 = get("BLOSM50_L1", envir = cachedEnv)$symbol,
        BLOSM50_L2 = get("BLOSM50_L2", envir = cachedEnv)$symbol,
        BLOSM50_L3 = get("BLOSM50_L3", envir = cachedEnv)$symbol,
        BLOSM50_L4 = get("BLOSM50_L4", envir = cachedEnv)$symbol,
        BLOSM50_L5 = get("BLOSM50_L5", envir = cachedEnv)$symbol,
        BLOSM50_L6 = get("BLOSM50_L6", envir = cachedEnv)$symbol,
        BLOSM50_L7 = get("BLOSM50_L7", envir = cachedEnv)$symbol,
        BLOSM50_L8 = get("BLOSM50_L8", envir = cachedEnv)$symbol,
        chemistry_property_Mahler = get("chemistry_property_Mahler", 
                                        envir = cachedEnv)$symbol, 
        contact_potential_Maiorov = get("contact_potential_Maiorov",
                                        envir = cachedEnv)$symbol, 
        protein_blocks_Rogov = get("protein_blocks_Rogov", 
                                   envir = cachedEnv)$symbol,
        sequence_alignment_Dayhoff = get("sequence_alignment_Dayhoff",
                                         envir = cachedEnv)$symbol, 
        structure_alignments_Mirny = get("structure_alignments_Mirny", 
                                         envir = cachedEnv)$symbol,
        custom = get("custom", envir = .globalEnv)$symbol)
}

#' Create sequence logo.
#' 
#' Create sequence logo for visualizing results of testing differential usage 
#' of amino acids.
#'
#' @param testDAUresults An object of \code{testDAUresults} (results of testing 
#' differential amino acid uasage).
#' @param type A character vector of length 1. Type of statistics to display 
#' on y-axis, available choices are "diff" or "zscore".
#' @param pvalueCutoff A numeric vector of length 1. A cutoff of p-value.
#' @param groupingSymbol A named character vector, with 3-letter symbols of amino acids
#' as names and single-letter symbols as values.
#' @param font A character vector of length 1. Font type for displaying sequence
#' Logo.
#' @param legend A logical vector of length 1, indicating whether to show the 
#' legend.
#' @param labelRelativeToAnchor A logical vector of length 1, indicating whether
#' x-axis label should be adjusted relative to the anchoring position.
#' @param labels A character vector, x-axis labels.
#' @param gscmd  A character vector of length 1, indicating the file path to 
#' Ghostscript. This should be defined as a system variable or provide directly.
#' @param fontface An integer, fontface of text for axis annotation and legends.
#' @param fontsize An integer, fontsize of text for axis annotation and legends.
#' 
#' @import grDevices
#' @import motifStack
#' @import grid
#' @importFrom grImport PostScriptTrace readPicture picture 
#'
#' @return A sequence Logo is plotted without returned values.
#' @export
#' @author Jianhong Ou, Haibo Liu
#'
#' @examples 
#' data('seq.example') 
#' data('proteome.example')
#' bg <- buildBackgroundModel(seq.example, proteome=proteome.example, 
#'                            numSubsamples=10, testType = "ztest")
#' t0 <- testDAU(seq.example, bg)
#' t1 <- testDAU(dagPeptides = seq.example, dagBackground = bg, 
#'               groupingScheme = "classic")
#' t2 <- testDAU(dagPeptides = seq.example, dagBackground = bg, 
#'              groupingScheme = "charge")
#' t3 <- testDAU(dagPeptides = seq.example, dagBackground = bg, 
#'              groupingScheme = "chemistry")
#' t4 <- testDAU(dagPeptides = seq.example, dagBackground = bg, 
#'               groupingScheme = "hydrophobicity")
#' dagLogo(t0)
#' dagLogo(t1, groupingSymbol = getGroupingSymbol(t1@group))
#' dagLogo(t2, groupingSymbol = getGroupingSymbol(t2@group))
#' dagLogo(t3, groupingSymbol = getGroupingSymbol(t3@group))
#' dagLogo(t4, groupingSymbol = getGroupingSymbol(t4@group))


dagLogo <- function(testDAUresults,
                    type = c("diff", "zscore"),
                    pvalueCutoff = 0.05,
                    groupingSymbol = getGroupingSymbol(testDAUresults@group),
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
    data <- getData(type, testDAUresults)
    dat <- data$dat
    npos <- ncol(dat)
    ncha <- nrow(dat)
    if (!is.null(labels)) 
    {
        if (length(labels) < npos) 
        {
            stop("The length of labels is too short!", call. = FALSE)
        }
    }
    colset <- colorsets(testDAUresults@group)
    colset <- colset[rownames(dat)]
    if (any(is.na(colset)))
        stop("Not every symbol has its color setting.", call. = FALSE)
    if (!is.null(groupingSymbol)) 
    {
        rownames(dat) <- groupingSymbol[rownames(dat)]
        names(colset) <- groupingSymbol[names(colset)]
    }
    rname <- rownames(dat)
    if (max(nchar(rname)) > 1) 
    {
        stop("Please using the groupingSymbol to convert the symbols into single letters.",
             call. = FALSE)
    }
    key <- paste("x",
                 ncha,
                 font,
                 paste(colset[rname], collapse = ""),
                 paste(rname, collapse = ""),
                 sep = "_")
    
    ## can't use ifelse here
    if(exists("tmp_motifStack_symbolsCache", envir = cachedEnv))
    {
        symbolsCache = get("tmp_motifStack_symbolsCache", envir = cachedEnv)   
    } else
    {
        symbolsCache = list()
    }
    if (!is.null(symbolsCache[[key]])) 
    {
        symbols <- symbolsCache[[key]]
    } else 
    {
        symbols <- motifStack:::coloredSymbols(ncha, font, colset[rname], rname)
        symbolsCache[[key]] <- symbols
        
        ## save symbolsCache to the environment variable for future use
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
    grid.text(label= data$label, 
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
        if (is.null(groupingSymbol)) 
        {
            groupingSymbol <- get("no", envir = cachedEnv)$symbol
        }
        for (i in 1:length(groupingSymbol)) 
        {
            grid.text(
                groupingSymbol[i],
                x = (npos+2)*dw,
                y = .95 - i * 0.05 * lwd/2,
                just = c(.5, .5),
                gp = gpar(col = colset[groupingSymbol[i]], fontsize=fontsize, fontface = fontface))
            grid.text(
                names(groupingSymbol)[i],
                x = (npos+2)*dw + 0.03 * lwd/2,
                y = .95 - i * 0.05 * lwd/2,
                just = c(0, .5),
                gp = gpar(col = colset[groupingSymbol[i]], fontsize=fontsize, fontface = fontface))
        }
    }
}