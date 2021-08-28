test_that("getting color palette works", {
  expect_length(.get_pal("libd_layer_colors"), 8)
  expect_length(.get_pal("Okabe-Ito"), 8)
  expect_length(.get_pal("navy"), 2)
  expect_length(.get_pal(c("gray90", "navy")), 2)
})
