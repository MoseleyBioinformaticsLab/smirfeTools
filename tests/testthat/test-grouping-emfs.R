context("test-grouping-emfs")

peak_list_2_df = function(peak_list){
  purrr::map2_df(peak_list, names(peak_list), function(.x, .y){
    data.frame(peak = .x, complete_EMF = .y, stringsAsFactors = FALSE)
  })
}

basic_peaks = list(X111 = paste0("X", seq(1, 3)),
                   X112 = paste0("X", seq(4, 6)),
                   X113 = paste0("X", seq(7, 10)))

test_that("no overlaps just returns the same thing", {
  basic_df = peak_list_2_df(basic_peaks)

  out_df = group_emfs_by_peaks(basic_df)

  expect_equal(length(unique(out_df$grouped_EMF)), 3)
  expect_equal(length(unique(out_df$complete_EMF)), 3)
  expect_equal(length(unique(out_df$peak)), 10)
})

# the most basic case, where we have a single completely overlapping case
doubled_peaks = list(X111 = paste0("X", seq(1, 3)),
             X112 = paste0("X", seq(1, 3)),
             X113 = paste0("X", seq(4, 10)))

test_that("collapsing true doubles works", {
  doubled_df = peak_list_2_df(doubled_peaks)

  out_df = group_emfs_by_peaks(doubled_df)

  expect_equal(length(unique(out_df$grouped_EMF)), 2)
  expect_equal(length(unique(out_df$complete_EMF)), 3)
  expect_equal(length(unique(out_df$peak)), 10)
})

# Now lets introduce a single peak that maps to two EMFs
mapped_two = list(X111 = paste0("X", seq(1, 3)),
                  X112 = paste0("X", c(1, seq(4, 6))),
                  X113 = paste0("X", seq(7, 10)))

test_that("collapsing with a single overlapped peak", {
  two_df = peak_list_2_df(mapped_two)

  out_df = group_emfs_by_peaks(two_df)

  expect_equal(length(unique(out_df$grouped_EMF)), 2)
  expect_equal(length(unique(out_df$complete_EMF)), 2)
  expect_equal(length(unique(out_df$peak)), 8)
})
