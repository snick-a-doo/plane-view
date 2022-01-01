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
    data <- c('start')
    for (i in 1:length(xss))
        data <- c(data, paste(xss[[i]]), '', paste(yss[[i]]), '')
    system2(pview.app, input=c(head(data, -1), 'end'))
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
    if (is.numeric(xs)) {
        if (nargs() == 1) { # 1. Y-vector given. Generate x.
            pvplot(list(1:length(xs)), list(xs))
            x.strs <- paste(1, length(xs), sep=':')
            y.strs <- as.character(call$xs)
        }
        else { # 2. X-vector and 1 or more y-vector given.
            yss <- list(...)
            pvplot(rep(list(xs), length(yss)), yss)
            x.strs <- as.character(call$xs)
            y.strs <- c$`...`
        }
    }
    else if (typeof(xs) == 'list') {
        if (nargs() == 1) { # 3. List of y-vectors given. Generate x for each.
            xss <- list()
            for (ys in xs)
                xss <- c(xss, list(1:length(xs[[1]])))
            pvplot(xss, xs)
        }
        else # 4. Lists of x-vectors and y-vectors given.
            pvplot(xs, list(...)[[1]])
    }
    else
        stop("Usage: Expecting 1 or more vectors or 1 or 2 lists.")

    ## !!Return a ggplot command that reproduces the plot.
    paste('ggplot() + geom_points(aes(', c$xs, ', ', c$`...`[[1]], '))')

}
