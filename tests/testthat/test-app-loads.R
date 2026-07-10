test_that("app.R defines a valid Shiny UI and server", {
  expect_true(inherits(app_env$ui, "shiny.tag") || inherits(app_env$ui, "shiny.tag.list"))
  expect_true(is.function(app_env$server))
  expect_equal(names(formals(app_env$server)), c("input", "output", "session"))
})

test_that("app.R's ui and server construct a valid Shiny app object", {
  app <- shiny::shinyApp(ui = app_env$ui, server = app_env$server)
  expect_s3_class(app, "shiny.appobj")
})
