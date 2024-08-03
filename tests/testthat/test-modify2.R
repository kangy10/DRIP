print("This is the beginning of test-modify2")

test_that("only accept proper input", {
  expect_error(modify2(bandwidth = 2, edge = 1:4),
               "edge must be a square matrix")
  expect_error(modify2(bandwidth = 2, edge = matrix(1.2, 4, 4), plot = FALSE),
               "edge can only have entry equal to 0 or 1")
  expect_error(modify2(bandwidth = -1, edge = matrix(1, 4, 4)),
               "bandwidth must be a positive integer")
  edge <- matrix(0, 9, 9)
  edge[4:5, ] <- 1
  edge[8, 8] <- 1
  expect_no_error(modify2(bandwidth = 1, edge = edge, plot = FALSE))
  edge <- matrix(0, 9, 9)
  edge[4:5, ] <- 1
  edge[8, 8] <- 1
  expect_no_error(modify2(bandwidth = 3, edge = edge, plot = TRUE))
  edge <- matrix(0, 9, 9)
  edge[1, 1] <- 1
  edge[9, 9] <- 1
  edge[1, 9] <- 1
  edge[9, 1] <- 1
  edge[5, 5] <- 1
  expect_no_error(modify2(bandwidth = 1, edge = edge, plot = TRUE))
  edge <- matrix(0, 9, 9)
  edge[1, 1] <- 1
  edge[9, 9] <- 1
  edge[1, 9] <- 1
  edge[9, 1] <- 1
  edge[5, 5] <- 1
  expect_no_error(modify2(bandwidth = 3, edge = edge, plot = TRUE))
  edge <- matrix(0, 9, 9)
  edge[1, 5] <- 1
  edge[9, 5] <- 1
  edge[5, 9] <- 1
  edge[5, 1] <- 1
  edge[5, 5] <- 1
  expect_no_error(modify2(bandwidth = 2, edge = edge, plot = TRUE))
  edge <- matrix(0, 9, 9)
  edge[1, 5] <- 1
  edge[9, 5] <- 1
  edge[5, 9] <- 1
  edge[5, 1] <- 1
  edge[5, 5] <- 1
  expect_no_error(modify2(bandwidth = 30, edge = edge, plot = TRUE))
  # edge <- matrix(0, 20, 20)
  # edge[2, 5] <- 1
  # edge[20, 5] <- 1
  # edge[5, 20] <- 1
  # edge[5, 2] <- 1
  # edge[5, 5] <- 1
  # expect_no_error(modify2(bandwidth = 2, edge = edge, plot = TRUE))
})

print("This is the end of test-modify2")
