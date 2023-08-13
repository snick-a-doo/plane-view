#' The name of the plain-view executable. Currently it's the full path to the development
#' build. Change to just the executable name after it's installed.
pview.app <- 'plane-view'

#' Helper function for pview()
#'
#' pvplot() spawns the plotting process and sends the data to it. It's not exported
#' because it requires a particular format for the data. pview() is a more flexible front
#' end to this function.
#'
#' @param xs An x-axis vector
#' @param yss A list of y-axis vectors
#'
#' @importFrom stringr str_split
#'
pvplot <- function(xs, yss, palette) {
    data <- c('pv.data')
    for (i in seq_along(yss))
        data <- c(data, paste(xs), 'pv.data', paste(yss[[i]]), 'pv.data')
    system2(pview.app,
            args = paste('--palette', palette),
            input = c(head(data, -1), 'pv.end'),
            stdout = TRUE) |>
        stringr::str_split(' ')
}

#' Interactive plotting
#'
#' @return
#' A string that gives the ranges of the plot.
#'
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom stringr str_glue
#'
#' @examples
#' wave <- tibble(theta = 0:99/10, cos = cos(theta), sin = sin(theta))
#' wave |> pview(theta, cos, sin)
#'
#' @export
pview <- function(.data, xs, ..., .palette = 'ggplot2') {
    xs <- .data |> dplyr::pull({{ xs }})
    yss <- purrr::map(enquos(...), function(ys) .data |> dplyr::pull({{ ys }}))
    range <- pvplot(xs, yss, .palette)[[1]]
    stringr::str_glue('coord_cartesian(xlim = c({range[1]}, {range[2]}), \\
                      ylim = c({range[3]}, {range[4]}))')
}
