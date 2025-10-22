test_that("get_chambers() returns numeric vector of chambers from an AquaResp experimental directory, not working now bc no ex data in pkg", {
  expect_equal(get_chambers(""), numeric(0))
})
