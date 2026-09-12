test_that("cluster_colours returns one colour per individual", {
  groups <- factor(c("A", "A", "B", "C", "B"))
  cols <- cluster_colours(groups)

  expect_length(cols, length(groups))
  expect_type(cols, "character")
})

test_that("individuals in the same group get the same colour", {
  groups <- factor(c("A", "B", "A", "B"))
  cols <- cluster_colours(groups)

  expect_identical(cols[1], cols[3])
  expect_identical(cols[2], cols[4])
  expect_false(cols[1] == cols[2])
})

test_that("colour assignment survives reordering", {
  groups <- factor(c("A", "A", "B", "C", "B"))
  cols <- cluster_colours(groups)

  idx <- c(5, 1, 3, 2, 4)
  reordered <- cluster_colours(groups[idx])

  # This is the whole point of the function: a positional rep() vector would
  # silently mislabel here, a lookup by group cannot.
  expect_identical(reordered, cols[idx])
})

test_that("a named palette is respected", {
  groups <- factor(c("A", "B", "A"))
  cols <- cluster_colours(groups, palette = c(A = "red", B = "blue"))

  expect_identical(unname(cols), c("red", "blue", "red"))
})

test_that("a named palette missing a group is an error", {
  groups <- factor(c("A", "B"))
  expect_error(
    cluster_colours(groups, palette = c(A = "red")),
    "No colour given"
  )
})

test_that("a short unnamed palette recycles with a warning", {
  groups <- factor(c("A", "B", "C"))
  expect_warning(
    cluster_colours(groups, palette = c("red", "blue")),
    "recycling"
  )
})

test_that("the legend attribute maps levels to colours", {
  groups <- factor(c("A", "B", "A"))
  cols <- cluster_colours(groups)
  legend <- attr(cols, "legend")

  expect_named(legend, c("A", "B"))
  expect_length(legend, 2)
})
