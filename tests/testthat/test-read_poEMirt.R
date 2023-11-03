# Load data
data(sim_data_dynamic)
data(sim_data_dynamic_repeated)

test_that(
  "Read data for static estimation",
  {
    data <- read_poEMirt(
      dataframe = sim_data_dynamic,
      response = paste0("y", 1:5),
      i = "i",
      j = "j"
    )
    expect_true(!exists("dynamic", data))
  }
)

test_that(
  "Read data for dynamic estimation",
  {
    data <- read_poEMirt(
      dataframe = sim_data_dynamic,
      response = paste0("y", 1:5),
      i = "i",
      j = "j",
      t = "t"
    )
    expect_true(exists("dynamic", data))
  }
)

test_that(
  "Read data with repeated items",
  {
    data <- read_poEMirt(
      dataframe = sim_data_dynamic_repeated,
      response = paste0("y", 1:5),
      i = "i",
      j = "j",
      t = "t"
    )
    expect_true(exists("rep", data))
    expect_identical(ncol(data$rep$raw), data$size$J)
    expect_identical(
      length(data$rep$processed), 
      length(unique(sim_data_dynamic_repeated$j))
    )
    expect_identical(
      length(unlist(data$rep$processed)),
      data$size$J
    )
  }
)