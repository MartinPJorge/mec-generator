#' @export
#' @title addLinkProps
#' @name addLinkProps
#' @description Includes additional link properties to the existing links
#' @param links data.frame of the links generated
#' @param from_ vector of ids of link origin within links$from
#' @param to_ vector of ids of link destination within links$to
#' @param properties named list with properties to be added
#' @return data.frame with the new links and properties
addLinkProps <- function(links, from_, to_, properties) {
  linksL <- as.list(links)

  # Seek the affected links
  affected <- c()
  for (row in 1:length(from_)) {
    currAffect <- which(links$from == from_[row] & links$to == to_[row])
    if (length(currAffect) > 0) {
      affected <- c(affected, currAffect)
    }
  }

  for (prop in attributes(properties)$name) {
    # Set the affected links to the stablished property
    property <- c()
    if (is.numeric(properties[[prop]])) {
      property <- rep(x = 0, times = nrow(links))
    } else {
      property <- rep(x = "", times = nrow(links))
    }

    # Property already present in the links data.frame
    if (prop %in% names(links)) {
      property <- as.vector(linksL[[prop]])
    }
    property[affected] <- properties[[prop]]
    linksL[[prop]] <- property
  }

  return(as.data.frame(linksL))
}


#' @export
#' @title addNodeProps
#' @name addNodeProps
#' @description Includes additional node properties to the existing nodes
#' @param nodes data.frame with the generated nodes
#' @param id_ node id within the data.frame
#' @param properties named list with properties to be added
#' @return data.frame with the new nodes and properties
addNodeProps <- function(nodes, id_, properties) {
  nodesL <- as.list(nodes)

  # Seek the affected nodes
  affected <- c()
  for (row in 1:length(id_)) {
    currAffect <- which(nodes$id == id_[row])
    if (length(currAffect) > 0) {
      affected <- c(affected, currAffect)
    }
  }

  for (prop in attributes(properties)$name) {
    # Set the affected nodes to the stablished property
    property <- c()
    if (is.numeric(properties[[prop]])) {
      property <- rep(x = 0, times = nrow(nodes))
    } else {
      property <- rep(x = "", times = nrow(nodes))
    }

    # Property already present in the nodes data.frame
    if (prop %in% names(nodes)) {
      property <- as.vector(nodesL[[prop]])
    }
    property[affected] <- properties[[prop]]
    nodesL[[prop]] <- property
  }

  return(as.data.frame(nodesL))
}
