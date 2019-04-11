context("Addition of extra properties")
library(mecgen)

test_that("Link properties addition works", {

  links <- data.frame(from = c('a', 'a2', 'a3'), to = c('b', 'b2', 'b3'),
                      bw = c(1, 2, 3))

  links2 <- addLinkProps(links = links, from_ = c('a', 'a3'),
                         to_ = c('b', 'b3'),
                         properties = list(bw = c(2, 4),
                                           radios = c("LTE", "GSM")))

  expect_equal(as.vector(links2$bw), c(2, 2, 4))
  expect_equal(as.vector(links2$radios), c("LTE", "", "GSM"))
})


test_that("Node properties addition works", {

  nodes <- data.frame(id = c('a', 'a2', 'a3'), bw = c(1, 2, 3))

  nodes2 <- addNodeProps(nodes = nodes, id_ = c('a', 'a3'),
                         properties = list(bw = c(2, 4),
                                           radios = c("LTE", "GSM")))

  expect_equal(as.vector(nodes2$bw), c(2, 2, 4))
  expect_equal(as.vector(nodes2$radios), c("LTE", "", "GSM"))
})
