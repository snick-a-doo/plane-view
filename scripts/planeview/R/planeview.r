#' The name of the plain-view executable. Currently it's the full path to the development
#' build. Change to just the executable name after it's installed.
pview.app <- '/home/samv/programs/plane-view/build/app/plane-view'

#' Helper function for pview()
#'
#' pvplot() spawns the plotting process and sends the data to it. It's not exported
#' because it requires a particular format for the data. pview() is a more flexible front
#' end to this function.
#'
#' @param xss A list of x-axis vectors
#' @param yss A list of y-axis vectors
pvplot <- function(xss, yss) {
    data <- c('pv.start')
    for (i in 1:length(xss))
        data <- c(data, paste(xss[[i]]), 'pv.sep', paste(yss[[i]]), 'pv.sep')
    strsplit(system2(pview.app, input=c(head(data, -1), 'pv.end'), stdout=TRUE), split=' ')[[1]]
}

#' Plot vectors
#'
#' pview may be called as
#' 1. pview(ys)  A vector of indices will be generated for x.
#' 2. pview(xs, ys1, ys2, ...)
#' 3. pview(list(ys1, ys2, ...))  A vector of indices will be generated for x.
#' 4. pview(list(xs1, xs2, ...), list(ys1, ys2, ...))
#'
#' @return
#' A string that evaluates to a ggplot of the displayed plot. \emph{Not implemented yet.}
#'
#' @examples
#' xs <- 0:99/10
#' pview(xs, cos(xs))
#' pview(xs, cos(xs), sin(xs))
#'
#' @export
pview <- function(xs, ...) {
    call <- match.call(expand.dots = FALSE)
    range <- character(0)
    if (is.numeric(xs)) {
        if (nargs() == 1) { # 1. Y-vector given. Generate x.
            range <- pvplot(list(1:length(xs)), list(xs))
            x.strs <- paste(1, length(xs), sep=':')
            y.strs <- deparse(call$xs)
        }
        else { # 2. X-vector and 1 or more y-vector given.
            yss <- list(...)
            range <- pvplot(rep(list(xs), length(yss)), yss)
            x.strs <- deparse(call$xs)
            y.strs <- call$`...`
        }
    }
    else if (typeof(xs) == 'list') {
        if (nargs() == 1) { # 3. List of y-vectors given. Generate x for each.
            xss <- list()
            for (ys in xs)
                xss <- c(xss, list(1:length(xs[[1]])))
            range <- pvplot(xss, xs)
        }
        else # 4. Lists of x-vectors and y-vectors given.
            range <- pvplot(xs, list(...)[[1]])
    }
    else
        stop("Usage: Expecting 1 or more vectors or 1 or 2 lists.")

    ## Return a ggplot command that reproduces the plot.
    ##! Not implemented for cases 3 or 4 yet.
    cols = as.list(1 + 1:length(y.strs))
    paste('ggplot()',
          paste('geom_point(aes(', x.strs, ', ', y.strs, '), color=', cols, ')', sep='', collapse=' + '),
          paste('coord_cartesian(xlim=c(', range[1], ', ', range[2],
                '), ylim=c(', range[3], ', ', range[4], '))', sep=''),
          "labs(x='', y='')",
          sep=' + ')
}
