app_env <- new.env(parent = globalenv())

app_lines <- readLines(testthat::test_path("..", "..", "app.R"))
launch_idx <- grep("^shinyApp\\(ui = ui, server = server\\)\\s*$", app_lines)
if (length(launch_idx) == 1) {
  app_lines <- app_lines[seq_len(launch_idx - 1)]
}
eval(parse(text = app_lines), envir = app_env)
