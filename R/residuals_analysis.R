residuals_analysis <-
  function(
    x,
    sampling_rate = 1,
    frequencies   = seq(0.01, 0.49, 0.01) * sampling_rate,
    nodes         = 2,
    keep_node     = nodes,
    filter_fun    = NULL,
    metric        = c("rmse", "sse", "mse"),
    return_ggplot = T
  ) {



    # load the necessary packages

    requireNamespace("ggplot2", "ggplot", "aes", "theme_classic", "labs",
                                "scale_colour_brewer", "geom_point",
                                "geom_line", "geom_text", "geom_vline")
    requireNamespace("signal", "butter", "filtfilt")
    requireNamespace("magrittr", "%>%")
    requireNamespace("purrr", "map_dbl", "pmap_dfr", "map_dfr")
    requireNamespace("dplyr", "mutate", "filter")
    requireNamespace("Rdpack", "reprompt")
    requireNamespace("stats", "lm", "predict", "setNames")
    requireNamespace("utils", "combn")



    # check the entered data

    stopifnot(is.numeric(x))
    stopifnot(is.numeric(sampling_rate))
    stopifnot(is.numeric(frequencies))
    stopifnot(is.numeric(nodes))
    stopifnot(is.numeric(keep_node))
    stopifnot(is.logical(return_ggplot))
    stopifnot(keep_node <= nodes)
    stopifnot(any(metric %in% c("rmse", "sse", "mse")))
    stopifnot(length(frequencies) >= nodes + 2)
    stopifnot(length(x) >= nodes + 2)

    if (is.null(filter_fun))
      filter_fun <-
        function(f, x) {
          flt <-
            butter(
              n     = 4,
              W     = f / (sampling_rate / 2),
              type  = "low",
              plane = "z"
            )
          return(filtfilt(flt, x))
        }

    nodes <-
      as.integer(nodes)

    keep_node <-
      as.integer(keep_node)

    if (length(metric) > 1)
      metric <-
        metric[1]


    if (metric == "mse")
      loss <-
        function(x, y) mean((x - y) ^ 2)

    else if (metric == "rmse")
      loss <-
        function(x, y) sqrt(mean((x - y) ^ 2))

     else if (metric == "sse")
      loss <-
        function(x, y) sum((x - y) ^ 2)



    # generate the map of residuals

    residuals <-
      map_dbl(frequencies, function(f) loss(x, filter_fun(f, x)))



    # calculate the SSE for each combination of nodes and sort them
    # NOTE: the edges of each segment overlaps

    fit_segment <-
      function(n1, n2) {

        idx <-
          seq(n1, n2)

        k <-
          frequencies[idx]

        y <-
          residuals[idx]

        z <-
          predict(lm(y ~ k))

        return(data.frame(x = k, y = y, z = z))
      }

    get_metric_df <-
      function(...) {

      df <-
        data.frame(...)

      metric <-
        map_dbl(1:(ncol(df) - 1), function(x) {

          f <-
            fit_segment(df[[names(df)[x]]], df[[names(df)[x + 1]]])

          if (x > 1) {
            f <-
              f[2:nrow(f),]
          }

          return(loss(f$y, f$z))
        }
        ) %>%
        sum()


      return(mutate(df, Metric = metric))
    }

    df <-
      seq(2, length(frequencies) - 1) %>%
      combn(., m = nodes) %>%
      t() %>%
      as.data.frame() %>%
      setNames(paste0("Node", seq(1, nodes))) %>%
      cbind(data.frame(Node0 = 1),
            .,
            data.frame(NodeN = length(frequencies))) %>%
      pmap_dfr(.l = ., .f = get_metric_df) %>%
      filter(Metric == min(Metric))



    # get the optimal frequency cutoff

    opt <-
      frequencies[df[[keep_node + 1]]]

    if (!return_ggplot)
      return(opt)



    # create the plot

    digits <-
      abs(min(0, floor(log10(max(frequencies))))) + 2

    label <-
      paste("Cut-off = ", format(opt, digits = digits, nsmall = digits), "Hz")

    opt.df <-
      data.frame(
        f      = opt,
        y      = 0.5 * (max(residuals) - min(residuals)),
        Legend = "Optimal frequency",
        l      = label
      )

    res.df <-
      data.frame(
        f = frequencies,
        y = residuals,
        Legend = "Residuals"
      )

    seg.df <-
      map_dfr(1:(nodes + 1), function(k) {

        f <-
          fit_segment(df[[names(df)[k]]], df[[names(df)[k + 1]]])

        return(
          data.frame(
            f = f$x,
            y = f$z,
            Legend = paste("Fitted segment ", k)
            )
          )
        }
      )

    pp <-
      ggplot(
        mapping = aes(
          x      = f,
          y      = y,
          colour = Legend
        )
      ) +
      theme_classic() +
      labs(
        x     = "Frequencies (Hz)",
        y     = metric,
        title = "Residuals Analysis"
      ) +
      scale_colour_brewer(palette = "Set1") +
      geom_point(
        data        = res.df,
        size        = 2,
        alpha       = 0.2
      ) +
      geom_line(
        data     = seg.df,
        size     = 1,
        linetype = "dashed"
      ) +
      geom_vline(
        data        = opt.df,
        mapping     = aes(xintercept = f),
        colour      = "black",
        linetype    = "solid",
        size        = 1.5,
        alpha       = 0.8
      ) +
      geom_text(
        data         = opt.df,
        mapping      = aes(label = l),
        hjust        = 0,
        size         = 4,
        colour       = "black",
        nudge_x      = 0.01 * (max(frequencies) - min(frequencies)),
        show.legend  = FALSE
      )

    return(list(opt = opt, plot = pp))
}
