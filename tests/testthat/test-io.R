test_that("read_genalex_csv fails clearly on a missing file", {
  expect_error(
    read_genalex_csv("definitely-not-here.csv"),
    "File not found"
  )
})

test_that("save_figure rejects an unknown extension", {
  expect_error(
    save_figure(tempfile(fileext = ".bmp"), plot(1:10)),
    "Unsupported file extension"
  )
})

test_that("save_figure writes a PDF", {
  f <- tempfile(fileext = ".pdf")
  save_figure(f, plot(1:10))

  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)

  unlink(f)
})
